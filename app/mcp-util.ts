/**
 * Helpers for the MCP server (app/mcp.ts): structured errors, value decoding,
 *  output summaries & step status. See docs/mcp.md for the contract these implement.
 */
import { z } from 'zod'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import type KRG from '@/core/KRG'
import type FPPRG from '@/core/FPPRG'
import type { FPL, Process } from '@/core/FPPRG'
import type { DataMetaNode, ProcessMetaNode, PromptMetaNode } from '@/spec/metanode'
import { decode_complete_process_inputs } from '@/core/engine'
import { TimeoutError } from '@/spec/error'

export const PUBLIC_URL = process.env.PUBLIC_URL || 'https://playbook-workflow-builder.cloud'
export const reportLink = (workflow_id: string) => `${PUBLIC_URL}/report/${workflow_id}`

/* ---------------------------------------------------------------- errors */

export type Issue = {
  /** where the problem is, e.g. `step.arguments.gene.step_id` */
  path: string
  /** machine readable problem, e.g. `type_mismatch` */
  code: string
  message: string
  expected?: unknown
  got?: unknown
  /** concrete values which would fix the problem */
  suggestions?: unknown[]
}

export class ToolError extends Error {
  constructor(public code: string, message: string, public issues: Issue[] = [], public hint?: string) {
    super(message)
  }
}

export const IssueC = z.object({
  path: z.string(),
  code: z.string(),
  message: z.string(),
  expected: z.unknown().optional(),
  got: z.unknown().optional(),
  suggestions: z.array(z.unknown()).optional(),
})
export const ErrorC = z.object({
  code: z.string(),
  message: z.string(),
  issues: z.array(IssueC).optional(),
  hint: z.string().optional(),
})
export type ToolErrorJSON = z.infer<typeof ErrorC>

export function toErrorJSON(exc: unknown): ToolErrorJSON {
  if (exc instanceof ToolError) {
    return { code: exc.code, message: exc.message, issues: exc.issues.length ? exc.issues : undefined, hint: exc.hint }
  }
  if (TimeoutError.isinstance(exc)) {
    return { code: 'timeout', message: 'Timed out waiting on the server, try again', hint: 'Retrying the same call is safe.' }
  }
  return { code: 'internal', message: exc instanceof Error ? exc.message : String(exc) }
}

/** Wrap a tool handler, reporting results and errors consistently */
export function tool<P, R extends Record<string, unknown>>(name: string, handler: (props: P) => Promise<R>) {
  return async (props: P) => {
    try {
      console.debug(JSON.stringify({ [name]: { input: props } }))
      const result = await handler(props)
      return {
        content: [{ type: 'text' as const, text: JSON.stringify(result) }],
        structuredContent: result,
      }
    } catch (exc) {
      const error = toErrorJSON(exc)
      if (error.code === 'internal') console.error(exc)
      console.error(JSON.stringify({ [name]: { input: props, error } }))
      return {
        content: [{ type: 'text' as const, text: JSON.stringify({ error }) }],
        structuredContent: { error },
        isError: true,
      }
    }
  }
}

/* ---------------------------------------------------------------- search */

function tokens(s: string) {
  return s.toLowerCase().split(/[^a-z0-9]+/).filter(t => t.length > 1)
}

/** A cheap relevance score: the number of query tokens found in the haystack */
export function score(q: string, haystack: (string | undefined)[]) {
  const hay = haystack.filter(Boolean).join(' ').toLowerCase()
  const hayTokens = new Set(tokens(hay))
  let n = 0
  for (const t of tokens(q)) {
    if (hayTokens.has(t)) n += 2
    else if (hay.includes(t)) n += 1
  }
  return n
}

export function processHaystack(proc: ProcessMetaNode) {
  return [
    proc.spec,
    proc.meta.label,
    proc.meta.description,
    proc.output.spec,
    proc.output.meta.label,
    ...dict.values(proc.inputs).map(i => array.ensureOne(i).spec),
    ...dict.values(proc.meta.tags ?? {}).flatMap(t => dict.keys(t) as string[]),
  ]
}

/** Find types close to an unknown one, for suggestions */
export function similarTypes(q: string, specs: { spec: string, meta: { label: string, description?: string } }[], n = 5) {
  return specs
    .map(node => ({ spec: node.spec, score: score(q, [node.spec, node.meta.label]) + (node.spec.toLowerCase().includes(q.toLowerCase()) ? 3 : 0) }))
    .filter(({ score }) => score > 0)
    .sort((a, b) => b.score - a.score)
    .slice(0, n)
    .map(({ spec }) => spec)
}

/* ---------------------------------------------------------------- summaries */

/** Make story text readable outside of the report: references become (doi:...), unbound placeholders are dropped */
export function cleanText(text: string | undefined) {
  if (!text) return ''
  const cleaned = text
    .replace(/\\ref\{([^}]*)\}/g, ' ($1)')
    .replace(/\\figref\{[^}]*\}/g, '')
    .replace(/\bundefined\b\s*/g, '')
    .replace(/\s+([.,;)])/g, '$1')
    .replace(/\s{2,}/g, ' ')
    .trim()
  return cleaned.charAt(0).toUpperCase() + cleaned.slice(1)
}

/**
 * Shrink a (decoded) value so that it fits into a model's context:
 *  long arrays, strings & objects are truncated with a marker.
 */
export function summarize(value: unknown, opts: { maxArray?: number, maxString?: number, maxKeys?: number, maxDepth?: number, marker?: boolean } = {}, depth = 0): unknown {
  const { maxArray = 10, maxString = 500, maxKeys = 30, maxDepth = 5, marker = true } = opts
  if (typeof value === 'string') {
    return value.length > maxString ? `${value.slice(0, maxString)}${marker ? `…(+${value.length - maxString} chars)` : ''}` : value
  }
  if (value === null || typeof value !== 'object') return value
  if (depth >= maxDepth) return Array.isArray(value) ? `[array of ${value.length}]` : '{…}'
  if (Array.isArray(value)) {
    const items = value.slice(0, maxArray).map(v => summarize(v, opts, depth + 1))
    if (marker && value.length > maxArray) items.push(`…(+${value.length - maxArray} more, ${value.length} total)`)
    return items
  }
  const entries = Object.entries(value as Record<string, unknown>)
  const out: Record<string, unknown> = {}
  for (const [k, v] of entries.slice(0, maxKeys)) out[k] = summarize(v, opts, depth + 1)
  if (marker && entries.length > maxKeys) out['…'] = `+${entries.length - maxKeys} more keys`
  return out
}

/** A short example which is still a valid value (no truncation markers) */
export function exampleOf(node: DataMetaNode | ProcessMetaNode | undefined) {
  if (!node) return undefined
  const example = node.meta.example ?? (node.kind === 'process' ? node.output.meta.example : undefined)
  if (example === undefined) return undefined
  return summarize(example, { maxArray: 5, maxString: 200, marker: false })
}

/** JSON Schema for a codec's zod spec where available */
export function jsonSchemaOf(codec: { zod?: z.ZodType } | undefined) {
  if (!codec?.zod) return undefined
  try {
    const { $schema: _, ...schema } = z.toJSONSchema(codec.zod, { unrepresentable: 'any', io: 'input' }) as Record<string, unknown>
    return schema
  } catch {
    return undefined
  }
}

export function isPrompt(node: ProcessMetaNode): node is PromptMetaNode {
  return 'prompt' in node && 'codec' in node
}

export function valueFormat(node: ProcessMetaNode) {
  if (!isPrompt(node)) return undefined
  const schema = jsonSchemaOf(node.codec)
  return {
    required: true,
    json_schema: schema ?? { type: 'string' },
    example: exampleOf(node),
  }
}

/* ---------------------------------------------------------------- values */

/**
 * Try to decode a value provided by the model with the prompt's codec, accepting
 *  some common alternative forms (JSON strings, bare lists for sets, ...).
 *  Returns the value as it would be stored by the UI (`codec.encode(decoded)`).
 */
export function decodeValue(node: PromptMetaNode, raw: unknown): { ok: true, encoded: unknown, decoded: unknown } | { ok: false, issue: Issue } {
  const candidates: unknown[] = [raw]
  if (typeof raw === 'string') {
    try { candidates.push(JSON.parse(raw)) } catch {}
  } else {
    candidates.push(JSON.stringify(raw))
  }
  // sets: allow a bare list, or a delimited string of terms
  if (Array.isArray(raw)) candidates.push({ set: raw.map(String) })
  if (typeof raw === 'string') {
    const set = raw.split(/[\n,;\t]+/).map(s => s.trim()).filter(Boolean)
    if (set.length > 0) candidates.push({ set })
  }
  let firstError: unknown
  for (const candidate of candidates) {
    try {
      const decoded = node.codec.decode(candidate as any)
      return { ok: true, decoded, encoded: node.codec.encode(decoded) }
    } catch (e) {
      if (firstError === undefined) firstError = e
    }
  }
  return {
    ok: false,
    issue: {
      path: 'step.value',
      code: 'invalid_value',
      message: `Value is not a valid ${node.meta.label}: ${formatDecodeError(firstError)}`,
      expected: jsonSchemaOf(node.codec) ?? { type: 'string' },
      suggestions: exampleOf(node) !== undefined ? [exampleOf(node)] : undefined,
    },
  }
}

function formatDecodeError(e: unknown) {
  if (e instanceof z.ZodError) {
    return e.issues.slice(0, 5).map(i => `${i.path.length ? i.path.join('.') : '(root)'}: ${i.message}`).join('; ')
  }
  return e instanceof Error ? e.message : String(e)
}

/* ---------------------------------------------------------------- arguments */

/** An argument reference as the model may provide it */
export const ArgRefC = z.union([
  z.string().describe('a step_id'),
  z.array(z.string()).describe('a list of step_ids (for arguments that take multiple)'),
  z.object({ step_id: z.string() }).strip(),
  z.object({ step_ids: z.array(z.string()) }).strip(),
  z.object({ ref: z.string().describe('(plan only) the ref of an earlier step in the same plan') }).strip(),
  z.object({ refs: z.array(z.string()) }).strip(),
])
export type ArgRef = z.infer<typeof ArgRefC>

export function normalizeArgRef(arg: ArgRef): { ids: string[], list: boolean } {
  if (typeof arg === 'string') return { ids: [arg], list: false }
  if (Array.isArray(arg)) return { ids: arg, list: true }
  if ('step_id' in arg) return { ids: [arg.step_id], list: false }
  if ('step_ids' in arg) return { ids: arg.step_ids, list: true }
  if ('ref' in arg) return { ids: [arg.ref], list: false }
  return { ids: arg.refs, list: true }
}

/** Describe an input of a process */
export function describeInput(input: DataMetaNode | DataMetaNode[]) {
  const node = array.ensureOne(input)
  return {
    type: node.spec,
    data_type: node.spec,
    multiple: Array.isArray(input),
    description: node.meta.description ?? node.meta.label,
  }
}

/** Process inputs are keyed `arg` or, for arguments taking multiple, `arg:0`, `arg:1`, ... */
export function processArguments(process: Process) {
  const args: Record<string, { step_id: string } | { step_ids: string[] }> = {}
  for (const { key, value } of dict.sortedItems(process.inputs)) {
    const [arg, ...i] = (key as string).split(':')
    if (i.length > 0) {
      const cur = args[arg]
      args[arg] = { step_ids: [...(cur && 'step_ids' in cur ? cur.step_ids : []), value.id] }
    } else {
      args[arg] = { step_id: value.id }
    }
  }
  return args
}

/* ---------------------------------------------------------------- status */

export type StepStatus =
  | { status: 'done', output?: unknown }
  | { status: 'error', error: string }
  | { status: 'waiting_for_input' }
  | { status: 'pending' }

/**
 * Get a step's status (and summarized output), waiting up to `waitMs` for it to resolve.
 *  Pending steps are always queued for resolution.
 */
export async function stepStatus(props: { krg: KRG, fpprg: FPPRG, process: Process, waitMs?: number, includeOutput?: boolean }): Promise<StepStatus> {
  const { krg, fpprg, process, waitMs = 0, includeOutput = false } = props
  const deadline = Date.now() + Math.min(waitMs, 60_000)
  let resolved = await fpprg.getResolved(process.id)
  if (!resolved) {
    if (waitMs <= 0) {
      // queue it without waiting
      fpprg.awaitResolved(process.id).catch(() => {})
    } else {
      while (!resolved && Date.now() < deadline) {
        try {
          resolved = await fpprg.awaitResolved(process.id)
        } catch (e) {
          if (!TimeoutError.isinstance(e)) throw e
        }
      }
    }
  }
  if (!resolved) return { status: 'pending' }
  if (!resolved.data) return { status: 'waiting_for_input' }
  if (resolved.data.type === 'Error') return { status: 'error', error: String(resolved.data.value) }
  if (!includeOutput) return { status: 'done' }
  let output: unknown = resolved.data.value
  try {
    const dataNode = krg.getDataNode(resolved.data.type)
    if (dataNode) output = dataNode.codec.decode(output as any)
  } catch {}
  return { status: 'done', output: summarize(output) }
}

/** A step's story abstract, using inputs/outputs only when they're already resolved (never blocks) */
export async function stepDescription(props: { krg: KRG, fpprg: FPPRG, process: Process }) {
  const { krg, fpprg, process } = props
  const metanode = krg.getProcessNode(process.type)
  if (!metanode) return ''
  let inputs: Record<string, unknown> | undefined
  let output: unknown
  try {
    const inputsResolved = await Promise.all(dict.values(process.inputs).map(p => fpprg.getResolved(p.id)))
    if (inputsResolved.every(r => r?.data !== undefined)) inputs = await decode_complete_process_inputs(krg, process)
  } catch {}
  try {
    const resolved = await fpprg.getResolved(process.id)
    if (resolved?.data && resolved.data.type === metanode.output.spec) output = metanode.output.codec.decode(resolved.data.value as any)
  } catch {}
  let data: unknown = process.data?.value
  try { if (isPrompt(metanode) && data !== undefined) data = metanode.codec.decode(data as any) } catch {}
  for (const props of [{ data, inputs, output }, { data }, {}]) {
    try {
      const abstract = metanode.story(props as any).abstract
      if (abstract) return cleanText(abstract)
    } catch {}
  }
  return metanode.meta.description ?? ''
}

/** Resolve the list of steps in a workflow */
export async function loadWorkflow(fpprg: FPPRG, workflow_id: string | null | undefined, path = 'workflow_id') {
  if (!workflow_id || workflow_id === 'start') return { fpl: undefined, steps: [] as FPL[] }
  const fpl = await fpprg.getFPL(workflow_id)
  if (!fpl) {
    throw new ToolError('workflow_not_found', `Workflow "${workflow_id}" not found`, [{
      path, code: 'not_found', message: `"${workflow_id}" is not a known workflow_id`,
    }], 'Use the workflow_id returned by your most recent expand/plan/set_value call, or omit it to start a new workflow.')
  }
  return { fpl, steps: fpl.resolve() }
}
