/**
 * The Playbook Workflow Builder MCP server, served (statelessly) at /mcp.
 *
 * See docs/mcp.md for the tool contract.
 */
import { McpServer } from '@modelcontextprotocol/sdk/server/mcp.js';
import { WebStandardStreamableHTTPServerTransport } from '@modelcontextprotocol/sdk/server/webStandardStreamableHttp.js';
import { getKRG } from '@/app/server-krg';
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { z } from 'zod'
import dedent from 'ts-dedent';
import fpprg from '@/app/fpprg';
import { FPL, type Process } from '@/core/FPPRG';
import type { ProcessMetaNode } from "@/spec/metanode"
import type KRG from "@/core/KRG"
import pluralize from "pluralize"
import { createCaller } from '@/app/trpc/router'
import { createServerContext } from '@/app/trpc/context'
import {
  ArgRefC, type ArgRef, ErrorC, cleanText, type Issue, ToolError,
  decodeValue, describeInput, exampleOf, isPrompt, jsonSchemaOf, loadWorkflow, normalizeArgRef,
  processArguments, processHaystack, reportLink, score, similarTypes, stepDescription, stepStatus,
  summarize, tool, valueFormat,
} from '@/app/mcp-util'

const trpc = createCaller(createServerContext())

/* ------------------------------------------------------------------ context */

/** A step output available to be used as an argument (an existing step, or an earlier step of a plan) */
type Available = { id: string, spec: string, label: string }

class StepContext {
  /** by step_id (or plan ref) */
  byId = new Map<string, Available>()
  /** most recent first */
  order: Available[] = []
  constructor(public krg: KRG) {}
  static fromWorkflow(krg: KRG, steps: FPL[]) {
    const ctx = new StepContext(krg)
    for (const step of steps) {
      const node = krg.getProcessNode(step.process.type)
      ctx.add({ id: step.process.id, spec: node?.output.spec ?? '?', label: node?.meta.label ?? step.process.type })
    }
    return ctx
  }
  add(avail: Available) {
    this.byId.set(avail.id, avail)
    this.order = [avail, ...this.order.filter(a => a.id !== avail.id)]
  }
  /** step_ids producing a given data type, most recent first (optionally preferring a given step) */
  candidates(spec: string, prefer?: string) {
    const ids = this.order.filter(a => a.spec === spec).map(a => a.id)
    if (prefer && ids.includes(prefer)) return [prefer, ...ids.filter(id => id !== prefer)]
    return ids
  }
}

/* ------------------------------------------------------------------ options */

/**
 * Not hidden, and the node registered under its spec: when two components share a spec, the KRG's
 *  input/output indices keep both but only the last one can be looked up (and so used).
 */
function usable(krg: KRG, proc: ProcessMetaNode) {
  return proc.meta.hidden !== true && krg.getProcessNode(proc.spec) === proc
}

function storyAbstract(proc: ProcessMetaNode) {
  try { return cleanText(proc.story({}).abstract ?? proc.meta.description) } catch { return cleanText(proc.meta.description) }
}

function describeOption(proc: ProcessMetaNode, ctx: StepContext, head?: string) {
  const args = dict.items(proc.inputs).map(({ key, value }) => {
    const input = describeInput(value)
    const candidates = ctx.candidates(input.data_type, head)
    const disabled_because = input.multiple
      ? candidates.length < 2 ? `Multiple ${pluralize(array.ensureOne(value).meta.label)} required (this workflow has ${candidates.length})` : undefined
      : candidates.length < 1 ? `${array.ensureOne(value).meta.label} required (no step in this workflow produces ${input.data_type})` : undefined
    return { key: key as string, input, candidates, disabled_because }
  })
  const disabled_because = args.map(a => a.disabled_because).filter(Boolean).join('; ') || undefined
  return {
    type: proc.spec,
    step_type: proc.spec,
    label: proc.meta.label,
    description: storyAbstract(proc) ?? '',
    output_type: proc.output.spec,
    arguments: dict.init(args.map(({ key, input, candidates }) => ({ key, value: { ...input, candidates } }))),
    suggested_arguments: disabled_because ? undefined : dict.init(args.map(({ key, input, candidates }) => ({
      key, value: input.multiple ? { step_ids: candidates } : { step_id: candidates[0] },
    }))),
    requires_value: isPrompt(proc) || undefined,
    value_example: isPrompt(proc) ? exampleOf(proc) : undefined,
    disabled_because,
    pagerank: proc.meta.pagerank ?? 0,
  }
}
type Option = ReturnType<typeof describeOption>

function computeOptions(props: { krg: KRG, steps: FPL[], head?: FPL, q?: string | null, limit?: number | null, include_disabled?: boolean | null }) {
  const { krg, steps, head } = props
  const ctx = StepContext.fromWorkflow(krg, steps)
  const headNode = head ? krg.getProcessNode(head.process.type) : undefined
  let items = krg.getNextProcess(headNode ? headNode.output.spec : '').filter(proc => usable(krg, proc))
  let options = items.map(proc => ({ option: describeOption(proc, ctx, head?.process.id), relevance: props.q ? score(props.q, processHaystack(proc)) : 0 }))
  if (props.q) options = options.filter(({ relevance }) => relevance > 0)
  if (props.include_disabled === false) options = options.filter(({ option }) => !option.disabled_because)
  options.sort((a, b) =>
    (b.relevance - a.relevance)
    || (Number(!a.option.disabled_because) > Number(!b.option.disabled_because) ? -1 : Number(!a.option.disabled_because) < Number(!b.option.disabled_because) ? 1 : 0)
    || (b.option.pagerank - a.option.pagerank)
  )
  const total = options.length
  const limit = props.limit ?? 30
  return { total, options: options.slice(0, limit).map(({ option: { pagerank: _, ...option } }) => option) }
}

/** A compact version of options to include in expand responses */
function compactOptions(options: Omit<Option, 'pagerank'>[]) {
  return options.filter(o => !o.disabled_because).slice(0, 15).map(o => ({
    type: o.type,
    label: o.label,
    output_type: o.output_type,
    suggested_arguments: o.suggested_arguments,
    requires_value: o.requires_value,
  }))
}

/* ------------------------------------------------------------------ steps */

const StepC = z.object({
  type: z.string().describe('the step type (`type` from options/find_steps), e.g. "Input[Gene]" or "GTExTissueExpressionFromGene" -- NOT a data type like "Term[Gene]"'),
  arguments: z.record(z.string(), ArgRefC).default({}).describe('for each argument in options: the step_id of a prior step producing that data type, e.g. {"gene": {"step_id": "..."}}; for arguments with multiple: true use {"step_ids": [...]}. options returns suggested_arguments you can pass as is.'),
  value: z.unknown().optional().describe('required for steps with requires_value (Input[...] steps & other prompts): the value, as described by describe(type).value; e.g. "ACE2" for Input[Gene], ["ACE2","TP53"] or {"description":"...","set":["ACE2","TP53"]} for Input[Set[Gene]]'),
})
type StepInput = z.infer<typeof StepC>

type ValidatedStep = {
  processNode: ProcessMetaNode
  inputs: Record<string, string>
  data?: { type: string, value: any }
  warnings: string[]
}

function unknownTypeIssue(krg: KRG, type: string, path: string): Issue {
  const dataNode = krg.getDataNode(type)
  if (dataNode) {
    const producers = krg.getPrevProcess(type).filter(p => p.meta.hidden !== true).map(p => p.spec)
    return {
      path, code: 'data_type_not_step_type',
      message: `"${type}" is a data type, not a step type. Use a step which produces it.`,
      suggestions: producers.slice(0, 10),
    }
  }
  return {
    path, code: 'unknown_step_type',
    message: `"${type}" is not a known step type. Use a type from options or find_steps.`,
    suggestions: similarTypes(type, krg.getProcessNodes().filter(p => p.meta.hidden !== true)),
  }
}

function validateStep(ctx: StepContext, step: StepInput, path = 'step'): ValidatedStep | { issues: Issue[] } {
  const { krg } = ctx
  const issues: Issue[] = []
  const warnings: string[] = []
  const processNode = krg.getProcessNode(step.type)
  if (!processNode) return { issues: [unknownTypeIssue(krg, step.type, `${path}.type`)] }
  const inputs: Record<string, string> = {}
  for (const argName in step.arguments) {
    const expected = processNode.inputs[argName]
    if (expected === undefined) {
      issues.push({
        path: `${path}.arguments.${argName}`, code: 'unexpected_argument',
        message: `${processNode.spec} has no argument "${argName}"`,
        suggestions: dict.keys(processNode.inputs) as string[],
      })
      continue
    }
    const input = describeInput(expected)
    const { ids, list } = normalizeArgRef(step.arguments[argName] as ArgRef)
    const candidates = ctx.candidates(input.data_type)
    if (!input.multiple && ids.length !== 1) {
      issues.push({
        path: `${path}.arguments.${argName}`, code: 'expected_single',
        message: `"${argName}" takes a single step_id ({"step_id": "..."}), got ${ids.length}`,
        suggestions: candidates.slice(0, 5).map(step_id => ({ step_id })),
      })
      continue
    }
    if (input.multiple && ids.length < 1) {
      issues.push({
        path: `${path}.arguments.${argName}`, code: 'expected_multiple',
        message: `"${argName}" takes a list of step_ids ({"step_ids": [...]})`,
        suggestions: candidates.length ? [{ step_ids: candidates }] : undefined,
      })
      continue
    }
    ids.forEach((id, i) => {
      const argPath = `${path}.arguments.${argName}${list ? `[${i}]` : ''}`
      const avail = ctx.byId.get(id)
      if (!avail) {
        issues.push({
          path: argPath, code: 'invalid_step_id',
          message: `"${id}" is not a step in this workflow${candidates.length ? `; steps producing ${input.data_type}: ${candidates.join(', ')}` : `, and no step in it produces ${input.data_type}`}`,
          expected: input.data_type,
          suggestions: candidates.length ? candidates : krg.getPrevProcess(input.data_type).filter(p => p.meta.hidden !== true).slice(0, 5).map(p => `add a ${p.spec} step first`),
        })
      } else if (avail.spec !== '?' && avail.spec !== input.data_type) {
        issues.push({
          path: argPath, code: 'type_mismatch',
          message: `step "${id}" (${avail.label}) produces ${avail.spec} but "${argName}" needs ${input.data_type}`,
          expected: input.data_type, got: avail.spec,
          suggestions: candidates.length ? candidates : undefined,
        })
      }
      inputs[input.multiple ? `${argName}:${i}` : argName] = id
    })
  }
  for (const argName in processNode.inputs) {
    if (step.arguments[argName] !== undefined) continue
    const input = describeInput(processNode.inputs[argName])
    const candidates = ctx.candidates(input.data_type)
    issues.push({
      path: `${path}.arguments.${argName}`, code: 'missing_argument',
      message: `missing argument "${argName}" (${input.data_type}${input.multiple ? ', multiple' : ''})`,
      expected: input.data_type,
      suggestions: candidates.length
        ? [input.multiple ? { step_ids: candidates } : { step_id: candidates[0] }]
        : krg.getPrevProcess(input.data_type).filter(p => p.meta.hidden !== true).slice(0, 5).map(p => `add a ${p.spec} step first`),
    })
  }
  let data: ValidatedStep['data']
  if (isPrompt(processNode)) {
    if (step.value === undefined || step.value === null || step.value === '') {
      issues.push({
        path: `${path}.value`, code: 'value_required',
        message: `${processNode.spec} requires a value`,
        expected: jsonSchemaOf(processNode.codec) ?? { type: 'string' },
        suggestions: exampleOf(processNode) !== undefined ? [exampleOf(processNode)] : undefined,
      })
    } else {
      const decoded = decodeValue(processNode, step.value)
      if (decoded.ok) data = { type: processNode.spec, value: decoded.encoded }
      else issues.push({ ...decoded.issue, path: `${path}.value` })
    }
  } else if (step.value !== undefined && step.value !== null && step.value !== '') {
    warnings.push(`${processNode.spec} doesn't take a value, it was ignored`)
  }
  if (issues.length > 0) return { issues }
  return { processNode, inputs, data, warnings }
}

/**
 * Validate (and unless dry_run, add) a list of steps to a workflow, atomically.
 *  Steps can refer to earlier steps of the same list by their `ref`.
 */
async function buildSteps(props: { krg: KRG, workflow_id?: string | null, steps: (StepInput & { ref?: string | null })[], dry_run?: boolean | null, path: (i: number) => string }) {
  const { krg } = props
  const { fpl: base, steps: existing } = await loadWorkflow(fpprg, props.workflow_id)
  const ctx = StepContext.fromWorkflow(krg, existing)
  const issues: Issue[] = []
  const validated: { ref: string, step: ValidatedStep }[] = []
  props.steps.forEach((step, i) => {
    const ref = step.ref ?? String(i)
    if (props.steps.length > 1 && ctx.byId.has(ref)) {
      issues.push({ path: `${props.path(i)}.ref`, code: 'duplicate_ref', message: `ref "${ref}" is already used` })
    }
    const result = validateStep(ctx, step, props.path(i))
    if ('issues' in result) {
      issues.push(...result.issues)
      ctx.add({ id: ref, spec: krg.getProcessNode(step.type)?.output.spec ?? '?', label: step.type })
    } else {
      validated.push({ ref, step: result })
      ctx.add({ id: ref, spec: result.processNode.output.spec, label: result.processNode.meta.label })
    }
  })
  if (issues.length > 0) {
    throw new ToolError('invalid_step', `${issues.length} problem${issues.length > 1 ? 's' : ''} with the requested step${props.steps.length > 1 ? 's' : ''}, nothing was added`, issues,
      'Fix the issues (see suggestions) and retry. Use dry_run: true to check without adding anything.')
  }
  if (props.dry_run) return { base, validated, built: undefined, fpl: undefined }
  const refToId: Record<string, string> = {}
  let fpl = base
  const built: { ref: string, process: Process, step: ValidatedStep }[] = []
  for (const { ref, step } of validated) {
    const process = await fpprg.resolveProcess({
      type: step.processNode.spec,
      inputs: dict.init(dict.items(step.inputs).map(({ key, value }) => ({ key, value: { id: refToId[value] ?? value } }))),
      data: step.data,
    })
    refToId[ref] = process.id
    fpl = fpl ? fpl.extend(process) : new FPL(process)
    built.push({ ref, process, step })
  }
  fpl = await fpprg.upsertFPL(fpl as FPL)
  return { base, validated, built, fpl }
}

const waitMs = (props: { include_output?: boolean | null, wait_seconds?: number | null }) =>
  1000 * (props.wait_seconds ?? (props.include_output ? 10 : 0))

async function stepResult(krg: KRG, process: Process, props: { include_output?: boolean | null, wait_seconds?: number | null }) {
  const processNode = krg.getProcessNode(process.type)
  if (!processNode) {
    return { step_id: process.id, step_type: process.type, type: '?', data_type: '?', label: process.type, description: 'This step type is no longer available', status: 'error' as const, error: `Unknown step type ${process.type}` }
  }
  const status = await stepStatus({ krg, fpprg, process, waitMs: waitMs(props), includeOutput: props.include_output ?? false })
  return {
    step_id: process.id,
    step_type: process.type,
    type: processNode.output.spec,
    data_type: processNode.output.spec,
    label: processNode.meta.label,
    description: await stepDescription({ krg, fpprg, process }),
    ...status,
  }
}

/* ------------------------------------------------------------------ docs */

const guide = dedent`
  # Playbook Workflow Builder (PWB) MCP

  A workflow is an ordered list of steps. Each step has a **step type** (e.g. "Input[Gene]", "GTExTissueExpressionFromGene")
  and produces a value of a **data type** (e.g. "Term[Gene]", "Scored[Tissue]"). Step arguments take the step_id of a
  prior step that produced the required data type.

  ## Tools
  - search_published(search): find published workflows; their workflow_id works with view/options/expand/set_value.
  - options(workflow_id?, step_id?, q?, limit?): what can be added next, from the latest step (or a given step).
    Each option gives \`suggested_arguments\` (ready to pass to expand), \`requires_value\` & \`value_example\`, or \`disabled_because\`.
  - find_steps(q?, input_type?, output_type?): search all step types, e.g. "what produces Set[Gene]?" or "enrichment".
  - describe(type?, q?): details for a step or data type (value format/example, producers, consumers); no type lists data types.
  - expand(workflow_id?, step, dry_run?, include_output?): add one step. Returns the new workflow_id and next options.
  - plan(workflow_id?, steps, dry_run?): add several steps atomically; steps refer to earlier ones with {"ref": "..."}.
  - set_value(workflow_id, step_id, value): change the value of an input step (e.g. run a published workflow on your gene).
  - view(workflow_id, include_output?): list the steps of a workflow with status (and outputs).

  ## Rules
  - Every change returns a NEW workflow_id; always continue from the most recent one. Older ids stay valid:
    building on one starts a separate branch, and steps added since are not part of it.
  - step_ids are stable identifiers of steps within a workflow (and its extensions).
  - Steps with requires_value need \`value\`: a string for single terms ("ACE2"), a list or {"description","set"} for sets.
  - Steps resolve asynchronously; status is one of pending, done, error, waiting_for_input. Pass include_output: true
    (optionally wait_seconds, max 60) to wait for and get a summary of the output.
  - Errors are structured: {error: {code, message, issues: [{path, code, message, expected, suggestions}], hint}};
    suggestions are concrete values to retry with.
  - Show the user the report link: ${reportLink('{workflow_id}')}
`

const instructions = dedent`
  Playbook Workflow Builder - Construct bioinformatics workflows step by step

  ## Key Capabilities
  - find information about genes, diseases, drugs, or other biomedical concepts
  - operate on gene sets, drug sets, data matrices and more
  - enrichment analysis, signature search, and other bioinformatics utilities
  - generate plots and construct persistent reports with full provenance

  ## Usage Patterns
  1. Check 'search_published' for a published workflow that already does what the user wants. If one does, 'view' it,
     and use 'set_value' to run it on the user's inputs.
  2. Otherwise start with 'options' (no workflow_id) to see input steps, and 'expand' with one of them and a 'value'.
  3. Continue with 'options' on the returned workflow_id (or the 'next_options' included in expand's response), and 'expand'
     with the chosen type and its 'suggested_arguments'. Use 'find_steps'/'describe' to look for a specific capability,
     and 'plan' to add several steps at once.
  4. Always continue from the latest workflow_id. Give the user the report link: ${reportLink('{workflow_id}')}

  Call 'help' for the full contract.

  ## Example
  \`\`\`
  {"type": "function_call", "name": "options", "arguments": {"q": "gene"}}
  {"type": "function_call_output", "output": {"result":[{"type":"Input[Gene]","label":"Gene Input","requires_value":true,"value_example":"ACE2","arguments":{}, ...}], ...}}
  {"type": "function_call", "name": "expand", "arguments": {"step":{"type":"Input[Gene]","value":"ACE2"}}}
  {"type": "function_call_output", "output": {"result":{"workflow_id":"38b60ad6","step_id":"2b976512","type":"Term[Gene]","status":"done","next_options":[{"type":"GTExTissueExpressionFromGene","suggested_arguments":{"gene":{"step_id":"2b976512"}}, ...}], ...}}}
  {"type": "function_call", "name": "expand", "arguments": {"workflow_id":"38b60ad6","step":{"type":"GTExTissueExpressionFromGene","arguments":{"gene":{"step_id":"2b976512"}}},"include_output":true}}
  {"type": "function_call_output", "output": {"result":{"workflow_id":"811b99c8","step_id":"965ac0ed","type":"Scored[Tissue]","status":"done","output":{...}, ...}}}
  {"type": "agent_message", "content": "I've created the report for you at <${reportLink('811b99c8')}>."}
  \`\`\`

  ## Important Notes
  PWB has a network of possibilities. The next possible option is not always what the user wants as an end product,
  but it might be a means to get there: use find_steps with output_type to find what produces what you need.
`

/* ------------------------------------------------------------------ schemas */

const OptionC = z.object({
  type: z.string().describe('the step type, use as step.type in expand'),
  step_type: z.string(),
  label: z.string(),
  description: z.string(),
  output_type: z.string().describe('the data type this step produces'),
  arguments: z.record(z.string(), z.object({
    type: z.string().describe('the data type this argument requires'),
    data_type: z.string(),
    multiple: z.boolean(),
    description: z.string(),
    candidates: z.array(z.string()).describe('step_ids in the workflow producing this data type, most relevant first'),
  })),
  suggested_arguments: z.record(z.string(), z.unknown()).optional().describe('arguments ready to pass to expand'),
  requires_value: z.boolean().optional(),
  value_example: z.unknown().optional(),
  disabled_because: z.string().optional().describe('why this option cannot be used yet'),
})
const CompactOptionC = z.object({
  type: z.string(), label: z.string(), output_type: z.string(),
  suggested_arguments: z.record(z.string(), z.unknown()).optional(),
  requires_value: z.boolean().optional(),
})
const StepResultC = z.object({
  step_id: z.string(),
  step_type: z.string(),
  type: z.string().describe('the data type this step produces'),
  data_type: z.string(),
  label: z.string(),
  description: z.string(),
  status: z.enum(['pending', 'done', 'error', 'waiting_for_input']),
  output: z.unknown().optional().describe('a summary of the output (with include_output)'),
  error: z.string().optional(),
})
const StepTypeC = z.object({
  type: z.string(), label: z.string(), description: z.string(), output_type: z.string(),
  arguments: z.record(z.string(), z.object({ type: z.string(), data_type: z.string(), multiple: z.boolean(), description: z.string() })),
  requires_value: z.boolean().optional(),
})
const outputOptions = {
  include_output: z.boolean().nullish().describe('wait for the step to resolve (up to wait_seconds, default 10) and include a summary of its output'),
  wait_seconds: z.number().min(0).max(60).nullish().describe('how long to wait for resolution'),
}

function describeStepType(proc: ProcessMetaNode) {
  return {
    type: proc.spec,
    label: proc.meta.label,
    description: storyAbstract(proc) ?? '',
    output_type: proc.output.spec,
    arguments: dict.init(dict.items(proc.inputs).map(({ key, value }) => ({ key: key as string, value: describeInput(value) }))),
    requires_value: isPrompt(proc) || undefined,
  }
}

function unknownDataType(krg: KRG, type: string, path: string) {
  return new ToolError('unknown_type', `"${type}" is not a known data type`, [{
    path, code: 'unknown_data_type',
    message: krg.getProcessNode(type) ? `"${type}" is a step type; its output type is "${krg.getProcessNode(type).output.spec}"` : `"${type}" is not a known data type`,
    suggestions: krg.getProcessNode(type) ? [krg.getProcessNode(type).output.spec] : similarTypes(type, krg.getDataNodes()),
  }], 'describe() with no type lists the data types.')
}

/* ------------------------------------------------------------------ server */

/** Build the server; a new one is created for each (stateless) request */
function buildServer() {
  const server = new McpServer({
    name: 'PWBMCP',
    version: '2.0.0',
    title: 'Playbook Workflow Builder',
    websiteUrl: process.env.PUBLIC_URL,
  }, { instructions })

  server.registerTool('help', {
    title: 'Help',
    description: 'The full contract for using these tools: concepts, tools, rules, error format.',
    inputSchema: {},
    outputSchema: { guide: z.string() },
    annotations: { readOnlyHint: true },
  }, tool('help', async () => ({ guide })))

  server.registerTool('options', {
    title: 'Options',
    description: 'Step types which can be added next, from the given step (default: the latest step) of a workflow, or input steps to start a new workflow. Each option includes suggested_arguments ready for expand, or why it is disabled.',
    inputSchema: {
      workflow_id: z.string().nullish().describe('Most recent workflow id (omit to start a new workflow)'),
      step_id: z.string().nullish().describe('Step to expand from (default: the latest step)'),
      q: z.string().nullish().describe('Only options matching these keywords (e.g. "enrichment")'),
      limit: z.number().int().min(1).max(200).nullish().describe('Maximum number of options (default 30)'),
      include_disabled: z.boolean().nullish().describe('Include options that cannot be used yet (default true)'),
    },
    outputSchema: {
      result: OptionC.array().optional(),
      workflow_id: z.string().nullish(),
      step_id: z.string().nullish(),
      total: z.number().optional().describe('total matching options, before limit'),
      error: ErrorC.optional(),
    },
    annotations: { readOnlyHint: true },
  }, tool('options', async (props) => {
    const { steps } = await loadWorkflow(fpprg, props.workflow_id)
    let head: FPL | undefined = steps[steps.length - 1]
    if (props.step_id && props.step_id !== 'start') {
      head = steps.find(step => step.process.id === props.step_id)
      if (!head) {
        throw new ToolError('step_not_found', `Step "${props.step_id}" is not in this workflow`, [{
          path: 'step_id', code: 'invalid_step_id', message: `"${props.step_id}" is not a step in workflow "${props.workflow_id ?? ''}"`,
          suggestions: steps.map(step => step.process.id).reverse(),
        }], props.workflow_id ? 'Use a step_id from view, or omit step_id to expand from the latest step.' : 'step_id requires a workflow_id.')
      }
    }
    const krg = await getKRG()
    const { total, options } = computeOptions({ krg, steps, head, q: props.q, limit: props.limit, include_disabled: props.include_disabled })
    return { result: options, workflow_id: props.workflow_id ?? null, step_id: head?.process.id ?? null, total }
  }))

  server.registerTool('find_steps', {
    title: 'Find step types',
    description: 'Search all step types by keyword and/or by the data types they take or produce, e.g. output_type "Set[Gene]" for every way to get a gene set. Use options to see which apply to a workflow.',
    inputSchema: {
      q: z.string().nullish().describe('keywords, e.g. "enrichment" or "expression"'),
      input_type: z.string().nullish().describe('only steps taking this data type (e.g. "Term[Gene]")'),
      output_type: z.string().nullish().describe('only steps producing this data type (e.g. "Set[Gene]")'),
      limit: z.number().int().min(1).max(100).nullish().describe('default 20'),
    },
    outputSchema: { result: StepTypeC.array().optional(), total: z.number().optional(), error: ErrorC.optional() },
    annotations: { readOnlyHint: true },
  }, tool('find_steps', async (props) => {
    const krg = await getKRG()
    if (props.input_type && !krg.getDataNode(props.input_type)) throw unknownDataType(krg, props.input_type, 'input_type')
    if (props.output_type && !krg.getDataNode(props.output_type)) throw unknownDataType(krg, props.output_type, 'output_type')
    let procs = props.input_type ? krg.getNextProcess(props.input_type)
      : props.output_type ? krg.getPrevProcess(props.output_type)
      : krg.getProcessNodes()
    procs = procs.filter(p => usable(krg, p))
    if (props.output_type) procs = procs.filter(p => p.output.spec === props.output_type)
    if (props.input_type) procs = procs.filter(p => dict.values(p.inputs).some(i => array.ensureOne(i).spec === props.input_type))
    let scored = procs.map(proc => ({ proc, relevance: props.q ? score(props.q, processHaystack(proc)) : 0 }))
    if (props.q) scored = scored.filter(({ relevance }) => relevance > 0)
    scored.sort((a, b) => (b.relevance - a.relevance) || ((b.proc.meta.pagerank ?? 0) - (a.proc.meta.pagerank ?? 0)))
    return { result: scored.slice(0, props.limit ?? 20).map(({ proc }) => describeStepType(proc)), total: scored.length }
  }))

  server.registerTool('describe', {
    title: 'Describe a type',
    description: 'Details for a step type (arguments, output, value format & example) or a data type (example, which steps produce/consume it). Without a type, lists data types (optionally filtered by q).',
    inputSchema: {
      type: z.string().nullish().describe('a step type (e.g. "Input[Set[Gene]]") or data type (e.g. "Set[Gene]")'),
      q: z.string().nullish().describe('when listing data types, only those matching these keywords'),
    },
    outputSchema: { result: z.unknown().optional(), error: ErrorC.optional() },
    annotations: { readOnlyHint: true },
  }, tool('describe', async (props) => {
    const krg = await getKRG()
    if (!props.type) {
      let types = krg.getDataNodes().filter(d => d.meta.hidden !== true).map(d => ({ d, relevance: props.q ? score(props.q as string, [d.spec, d.meta.label, d.meta.description]) : 0 }))
      if (props.q) types = types.filter(({ relevance }) => relevance > 0).sort((a, b) => b.relevance - a.relevance)
      return { result: types.map(({ d }) => ({ type: d.spec, label: d.meta.label, description: cleanText(d.meta.description) })) }
    }
    const proc = krg.getProcessNode(props.type)
    if (proc) {
      return { result: { kind: 'step', ...describeStepType(proc), value: valueFormat(proc) } }
    }
    const data = krg.getDataNode(props.type)
    if (data) {
      const producers = krg.getPrevProcess(data.spec).filter(p => p.meta.hidden !== true)
      const consumers = krg.getNextProcess(data.spec).filter(p => p.meta.hidden !== true)
      return {
        result: {
          kind: 'data', type: data.spec, label: data.meta.label, description: data.meta.description,
          example: exampleOf(data), json_schema: jsonSchemaOf(data.codec),
          produced_by: producers.slice(0, 30).map(p => p.spec), n_produced_by: producers.length,
          consumed_by: consumers.slice(0, 30).map(p => p.spec), n_consumed_by: consumers.length,
        },
      }
    }
    throw new ToolError('unknown_type', `"${props.type}" is not a known step or data type`, [{
      path: 'type', code: 'unknown_type', message: `"${props.type}" is not a known step or data type`,
      suggestions: similarTypes(props.type, [...krg.getProcessNodes(), ...krg.getDataNodes()].filter(n => n.meta.hidden !== true)),
    }])
  }))

  server.registerTool('expand', {
    title: 'Create or Expand Workflow',
    description: 'Add a step to a workflow (or start a new one). Returns the NEW workflow_id (use it from now on), the step_id, its status and next options. Use dry_run to only validate.',
    inputSchema: {
      workflow_id: z.string().nullish().describe('Most recent workflow id (omit for a new workflow)'),
      step: StepC,
      dry_run: z.boolean().nullish().describe('only validate the step, add nothing'),
      ...outputOptions,
    },
    outputSchema: {
      result: StepResultC.extend({
        workflow_id: z.string().describe('the new workflow id, use it for every subsequent call'),
        parent_workflow_id: z.string().nullish(),
        link: z.string().describe('report link for the user'),
        warnings: z.array(z.string()).optional(),
        next_options: CompactOptionC.array().optional(),
      }).optional(),
      dry_run: z.object({ valid: z.boolean(), type: z.string(), data_type: z.string(), arguments: z.record(z.string(), z.string()) }).optional(),
      error: ErrorC.optional(),
    },
  }, tool('expand', async (props) => {
    const krg = await getKRG()
    if (!props.workflow_id && dict.keys(props.step.arguments ?? {}).length > 0) {
      throw new ToolError('workflow_required', 'Arguments refer to prior steps, which requires a workflow_id', [{
        path: 'workflow_id', code: 'missing', message: 'workflow_id is required when step.arguments are given',
      }], 'Pass the workflow_id returned by your most recent expand.')
    }
    const { base, validated, built, fpl } = await buildSteps({ krg, workflow_id: props.workflow_id, steps: [props.step], dry_run: props.dry_run, path: () => 'step' })
    if (!built || !fpl) {
      const step = validated[0].step
      return { dry_run: { valid: true, type: step.processNode.spec, data_type: step.processNode.output.spec, arguments: step.inputs } }
    }
    const { process, step } = built[0]
    const { options } = computeOptions({ krg, steps: fpl.resolve(), head: fpl })
    return {
      result: {
        workflow_id: fpl.id,
        parent_workflow_id: base?.id ?? null,
        link: reportLink(fpl.id),
        ...await stepResult(krg, process, props),
        warnings: step.warnings.length ? step.warnings : undefined,
        next_options: compactOptions(options),
      },
    }
  }))

  server.registerTool('plan', {
    title: 'Add several steps',
    description: 'Validate and add several steps at once (atomically: nothing is added if any step is invalid). Later steps refer to earlier ones with {"ref": "<ref>"} (ref defaults to the step\'s index, "0", "1", ...), or to existing steps with {"step_id": "..."}.',
    inputSchema: {
      workflow_id: z.string().nullish().describe('Most recent workflow id (omit for a new workflow)'),
      steps: z.array(StepC.extend({ ref: z.string().nullish().describe('a name for this step, for later steps to refer to') })).min(1).max(50),
      dry_run: z.boolean().nullish().describe('only validate the steps, add nothing'),
      ...outputOptions,
    },
    outputSchema: {
      result: z.object({
        workflow_id: z.string(),
        parent_workflow_id: z.string().nullish(),
        link: z.string(),
        steps: StepResultC.extend({ ref: z.string() }).array(),
        next_options: CompactOptionC.array().optional(),
      }).optional(),
      dry_run: z.object({ valid: z.boolean(), steps: z.array(z.object({ ref: z.string(), type: z.string(), data_type: z.string() })) }).optional(),
      error: ErrorC.optional(),
    },
  }, tool('plan', async (props) => {
    const krg = await getKRG()
    const { base, validated, built, fpl } = await buildSteps({ krg, workflow_id: props.workflow_id, steps: props.steps, dry_run: props.dry_run, path: i => `steps[${i}]` })
    if (!built || !fpl) {
      return { dry_run: { valid: true, steps: validated.map(({ ref, step }) => ({ ref, type: step.processNode.spec, data_type: step.processNode.output.spec })) } }
    }
    // only wait on the last step: its resolution resolves the others it depends on
    const steps = await Promise.all(built.map(async ({ ref, process }, i) => ({
      ref, ...await stepResult(krg, process, i === built.length - 1 ? props : { include_output: props.include_output, wait_seconds: 0 }),
    })))
    const { options } = computeOptions({ krg, steps: fpl.resolve(), head: fpl })
    return {
      result: {
        workflow_id: fpl.id,
        parent_workflow_id: base?.id ?? null,
        link: reportLink(fpl.id),
        steps,
        next_options: compactOptions(options),
      },
    }
  }))

  server.registerTool('set_value', {
    title: 'Change an input value',
    description: 'Change the value of a step with a value (e.g. the gene of an Input[Gene] step), recomputing the steps that depend on it. Use this to run a published workflow on different inputs. Returns the new workflow_id and the new step_ids.',
    inputSchema: {
      workflow_id: z.string().describe('the workflow to change'),
      step_id: z.string().describe('the step whose value to change (see view)'),
      value: z.unknown().describe('the new value, see describe(step_type)'),
      ...outputOptions,
    },
    outputSchema: {
      result: StepResultC.extend({
        workflow_id: z.string(),
        parent_workflow_id: z.string(),
        link: z.string(),
        step_id_changes: z.record(z.string(), z.string()).describe('old step_id -> new step_id, for every step that changed'),
      }).optional(),
      error: ErrorC.optional(),
    },
  }, tool('set_value', async (props) => {
    const krg = await getKRG()
    const { fpl, steps } = await loadWorkflow(fpprg, props.workflow_id)
    const old = steps.find(step => step.process.id === props.step_id)
    if (!fpl || !old) {
      throw new ToolError('step_not_found', `Step "${props.step_id}" is not in this workflow`, [{
        path: 'step_id', code: 'invalid_step_id', message: `"${props.step_id}" is not a step in workflow "${props.workflow_id}"`,
        suggestions: steps.filter(step => { const n = krg.getProcessNode(step.process.type); return n && isPrompt(n) }).map(step => step.process.id),
      }], 'Use view to see the steps (and which have a value).')
    }
    const processNode = krg.getProcessNode(old.process.type)
    if (!processNode || !isPrompt(processNode)) {
      throw new ToolError('no_value', `${old.process.type} doesn't take a value`, [{
        path: 'step_id', code: 'no_value', message: `${old.process.type} doesn't take a value`,
        suggestions: steps.filter(step => { const n = krg.getProcessNode(step.process.type); return n && isPrompt(n) }).map(step => step.process.id),
      }])
    }
    const decoded = decodeValue(processNode, props.value)
    if (!decoded.ok) throw new ToolError('invalid_value', decoded.issue.message, [{ ...decoded.issue, path: 'value' }])
    const newProcess = await fpprg.resolveProcess({
      type: old.process.type,
      inputs: dict.init(dict.items(old.process.inputs).map(({ key, value }) => ({ key: key as string, value: { id: value.id } }))),
      data: { type: processNode.spec, value: decoded.encoded as any },
    })
    const { head } = fpl.rebase(old.process, newProcess)
    const newFpl = await fpprg.upsertFPL(head)
    const newSteps = newFpl.resolve()
    const step_id_changes = dict.init(steps.flatMap((step, i) =>
      newSteps[i] && newSteps[i].process.id !== step.process.id ? [{ key: step.process.id, value: newSteps[i].process.id }] : []))
    return {
      result: {
        workflow_id: newFpl.id,
        parent_workflow_id: fpl.id,
        link: reportLink(newFpl.id),
        ...await stepResult(krg, newFpl.process, props),
        step_id_changes,
      },
    }
  }))

  server.registerTool('search_published', {
    title: 'Search Published Workflow Apps',
    description: 'Find published workflows that already do what the user wants. Their workflow_id works with view, options, expand and set_value (to run it on the user\'s inputs).',
    inputSchema: {
      search: z.string().describe('A search that will identify relevant workflows, structure of the workflow like tools, databases, input/output types are relevant, specific biological entities are NOT.'),
    },
    outputSchema: {
      result: z.object({
        workflow_id: z.string(), title: z.string(), description: z.string().nullish(),
        dataSources: z.string().nullish(), inputs: z.string().nullish(), outputs: z.string().nullish(), link: z.string(),
      }).array().optional(),
      error: ErrorC.optional(),
    },
    annotations: { readOnlyHint: true },
  }, tool('search_published', async (props) => {
    const [published, community] = await Promise.all([
      trpc.playbooks.public({ search: props.search, limit: 10 }),
      trpc.playbooks.community({ search: props.search, limit: 10 }),
    ])
    const result = [
      ...published.filter(p => !p.disabled).map(p => ({
        workflow_id: p.id, title: p.label, description: p.description, dataSources: p.dataSources, inputs: p.inputs, outputs: p.outputs,
      })),
      ...community.map(p => ({
        workflow_id: p.playbook, title: p.title ?? 'Untitled', description: p.description, dataSources: undefined, inputs: p.inputs, outputs: p.outputs,
      })),
    ].map(p => ({ ...p, link: reportLink(p.workflow_id) }))
    return { result }
  }))

  server.registerTool('view', {
    title: 'View a workflow by id',
    description: 'The steps of a workflow: step_ids, types, arguments, values, status, and (with include_output) output summaries.',
    inputSchema: {
      workflow_id: z.string(),
      step_id: z.string().nullish().describe('only include the output of this step'),
      ...outputOptions,
    },
    outputSchema: {
      result: z.object({
        workflow_id: z.string(),
        link: z.string(),
        title: z.string().nullish(),
        description: z.string().nullish(),
        steps: StepResultC.extend({
          workflow_id: z.string().describe('the workflow id ending at this step (to branch off from here)'),
          arguments: z.record(z.string(), z.unknown()),
          value: z.unknown().optional(),
        }).array(),
      }).optional(),
      error: ErrorC.optional(),
    },
    annotations: { readOnlyHint: true },
  }, tool('view', async (props) => {
    const { fpl, steps } = await loadWorkflow(fpprg, props.workflow_id)
    if (!fpl) throw new ToolError('workflow_not_found', 'workflow_id is required', [{ path: 'workflow_id', code: 'missing', message: 'workflow_id is required' }])
    if (props.step_id && !steps.some(step => step.process.id === props.step_id)) {
      throw new ToolError('step_not_found', `Step "${props.step_id}" is not in this workflow`, [{
        path: 'step_id', code: 'invalid_step_id', message: `"${props.step_id}" is not a step in this workflow`, suggestions: steps.map(s => s.process.id),
      }])
    }
    const krg = await getKRG()
    const result = await Promise.all(steps.map(async (step) => {
      const processNode = krg.getProcessNode(step.process.type)
      const outputProps = props.step_id && props.step_id !== step.process.id ? { include_output: false, wait_seconds: 0 } : props
      let value: unknown = step.process.data?.value
      if (value !== undefined) {
        try { if (processNode && isPrompt(processNode)) value = processNode.codec.decode(value as any) } catch {}
        value = summarize(value)
      }
      return {
        workflow_id: step.id,
        ...await stepResult(krg, step.process, outputProps),
        arguments: processArguments(step.process),
        value,
      }
    }))
    return {
      result: {
        workflow_id: fpl.id,
        link: reportLink(fpl.id),
        title: fpl.playbook_metadata?.title || null,
        description: fpl.playbook_metadata?.description || null,
        steps: result,
      },
    }
  }))

  return server
}

export async function handleMcp(req: Request) {
  // stateless: a server & transport per request (a server can only be connected to one transport at a time)
  const server = buildServer()
  const transport = new WebStandardStreamableHTTPServerTransport({
    sessionIdGenerator: undefined,
    enableJsonResponse: true,
  })
  await server.connect(transport)
  try {
    return await transport.handleRequest(req)
  } finally {
    await transport.close()
    await server.close()
  }
}
