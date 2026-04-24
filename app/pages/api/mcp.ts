import handler from '@/utils/next-rest'
import { McpServer } from '@modelcontextprotocol/sdk/server/mcp.js';
import { StreamableHTTPServerTransport } from '@modelcontextprotocol/sdk/server/streamableHttp.js';
import cache from '@/utils/global_cache';
import krg from '@/app/krg';
import * as dict from '@/utils/dict'
import { z } from 'zod'
import dedent from 'ts-dedent';
import fpprg from '@/app/fpprg';
import { FPL, Process } from '@/core/FPPRG';
import type { ProcessMetaNode } from "@/spec/metanode"
import type KRG from "@/core/KRG"
import pluralize from "pluralize"

function options(props: { krg: KRG, heads: FPL[], workflow: FPL[] }) {
  let items: ProcessMetaNode[]
  const processNode = props.heads[0] ? props.krg.getProcessNode(props.heads[0].process.type) : undefined
  // we'll use leaf nodes of the metapath + the current selected node as the selections
  const selections: Record<string, { process: Process, processNode: ProcessMetaNode }> = {}
  if (props.heads.length <= 1) {
    ;[...props.workflow, ...props.heads].forEach(item => {
      if (item === undefined) return
      // add this to the selections
      selections[item.process.id] = { process: item.process, processNode: props.krg.getProcessNode(item.process.type) }
      // if a selection previously registered is a parent of this selection, remove it from selections
      dict.values(item.process.inputs).forEach(k => {
        if (k.id in selections) delete selections[k.id]
      })
    })
    items = props.krg.getNextProcess(processNode ? processNode.output.spec : '')
  } else {
    props.heads.forEach(item => {
      if (item === undefined) return
      selections[item.process.id] = { process: item.process, processNode: props.krg.getProcessNode(item.process.type) }
    })
    items = props.krg.getNextProcess(processNode ? processNode.output.spec : '')
      .filter(proc => dict.values(proc.inputs).some(value => Array.isArray(value)))
  }
  items = items.filter(proc => proc.meta.hidden !== true)
  items.forEach(item => {
    // determine if multi-inputs are satisfiable, and record a reason if not
    dict.items(item.inputs).map(({ key, value }) => {
      if (Array.isArray(value)) {
        return {
          key,
          value: dict.values(selections).filter(selection => selection.processNode.output.spec === value[0].spec).length <= 1 ? `Multiple ${pluralize(value[0].meta.label)} required` : undefined,
        }
      } else {
        return {
          key,
          value: dict.values(selections).filter(selection => selection.processNode.output.spec === value.spec).length < 1 ? `${value.meta.label} required` : undefined,
        }
      }
    }).forEach(({ key, value }) => {
      if (value !== undefined) {
        Object.assign(item.meta, { disabledBecause: value })
      }
    })
  })
  return { selections, items }
}

const server = cache('mcp', () => {
  const server = new McpServer({
    name: 'PWBMCP',
    version: '1.0.0',
    title: 'Playbook Workflow Builder',
    websiteUrl: process.env.PUBLIC_URL,
  }, {
    instructions:  dedent`
      Playbook Workflow Builder - Construct bioinformatics workflows step by step

      ## Key Capabilities
      - find information about genes, diseases, or other biomedical concepts
      - operate on gene sets, drug sets, data matrices and more
      - facilitate enrichment analysis, signature search, and other bioinformatics utilities
      - generate plots and construct persistent reports with full provenance

      ## Usage Patterns

      At each step, use the 'options' tool with the current workflow and choose from the next possible options with confirmation from the user by using the 'expand' tool.
      To start a workflow, use the 'options' without specifying any 'workflow_id' or 'step_id', and subsequently 'expand' with the relevant step without a 'workflow_id' or 'step_id'.
      Step 'type' specified in 'expand' MUST appear in 'options'.
      The workflow_id can be used to build a link for the user with: https://playbook-workflow-builder.cloud/report/{workflow_id}

      ## Example

      {"type": "user_message", "content": "I'd like to get the expression of ACE2."}
      {"type": "function_call", "name": "options", "arguments": {}}
      {"type": "function_call_output", "output": {"result":[{"type":"Input[Gene]","derived_from":{},"label":"Gene Input","description":"The workflow starts with selecting a search term."}]}}
      {"type": "user_message", "content": "I'd like to get the expression of ACE2."}
      {"type": "function_call", "name": "expand", "arguments": {"step":{"type":"Input[Gene]","derived_from":{},"value":"ACE2"}}}
      {"type": "function_call_output", "output": {"result":{"workflow_id": "38b60ad6", "step_id": "2b976512", "type": "Gene"}}
      {"type": "function_call", "name": "options", "arguments": {"workflow_id":"38b60ad6", "step_id": "2b976512"}}
      {"type": "function_call_output", "output": {"result":[{"type":"GetTissueExpressionOfGene","derived_from":{"gene": {"type":"Gene"}},"label":"Get expression of the gene", ...}]}}
      {"type": "function_call", "name": "expand", "arguments": {"workflow_id":"38b60ad6", "step":{"type":"GetTissueExpressionOfGene","derived_from":{"gene": {"step_id": "2b976512"}}}}}
      {"type": "function_call_output", "output": {"result":{"workflow_id": "811b99c8", "step_id": "965ac0ed", "type": "TissueGeneExpression"}}
      {"type": "agent_message", "content": "I've created the report for you at https://playbook-workflow-builder.cloud/report/811b99c8 , let me know what else you'd like to do."}

      ## Important Notes

      PWB has a network of possibilities. The next possible option is not always what the user wants as an end product, but it might be a means to get there.
    `,
  })
  server.registerTool('options', {
    title: 'Options to start or expand a workflow with',
    inputSchema: {
      workflow_id: z.string().nullish().describe('The complete workflow we are building'),
      step_id: z.string().nullish().describe('A specific step in the workflow we wish to expand from'),
    },
    outputSchema: {
      result: z.object({ type: z.string(), label: z.string(), description: z.string(), derived_from: z.record(z.string(), z.object({ type: z.string(), description: z.string() }).describe('type of output this can be derived from')) }).array().optional(),
      error: z.any().or(z.undefined()),
    },
  }, async (props) => {
    try {
      if (props.workflow_id === 'start') props.workflow_id = undefined
      if (props.step_id === 'start') props.step_id = undefined
      const fpl = props.workflow_id ? await fpprg.getFPL(props.workflow_id) : undefined
      if (props.workflow_id && !fpl) throw new Error(`{"workflow_id":>"${props.workflow_id}" not found<}`)
      const workflow = fpl?.resolve() ?? []
      let head
      if (props.step_id) {
        for (const item of workflow ?? []) {
          if (item.process.id === props.step_id) {
            head = item
            break
          }
        }
        if (!head) throw new Error(`{"step_id":>"${props.workflow_id}" not found<}`)
      } else {
        head = workflow[workflow.length-1]
      }
      const { items } = options({ krg, heads: [head], workflow })
      const result = items.map(proc => ({
        type: proc.spec,
        derived_from: dict.init(dict.items(proc.inputs).map(({ key, value }) => ({ key, value: Array.isArray(value) ? { type: value[0].spec, description: value[0].meta.description } : { type: value.spec, description: value.meta.description } }))),
        label: proc.meta.label,
        description: proc.story({}).abstract ?? proc.meta.description,
      }))
      return {
        content: [{ type: 'text', text: JSON.stringify({ result }) }],
        structuredContent: { result },
      }
    } catch (exc) {
      const error = (exc as Error).message
      console.error(JSON.stringify({ options: { input: props, error } }))
      return {
        content: [{ type: 'text', text: JSON.stringify({ error }) }],
        structuredContent: { error },
      }
    }
  })
  server.registerTool('expand', {
    title: 'Expand Workflow',
    inputSchema: {
      workflow_id: z.string().nullish().describe('The complete workflow we are building off of (undefined when building a new workflow)'),
      step: z.object({
        type: z.string().describe('The type of operation as seen in `options`'),
        derived_from: z.record(z.string(), z.object({ step_id: z.string() }).strip()).default({}).describe('output of previous step_id to provide as input for this step'),
        value: z.string().nullish().describe('The input value for this step when not coming from a prior step'),
      }).strict().describe('What we do at this step'),
    },
    outputSchema: {
      result: z.object({ workflow_id: z.string(), step_id: z.string(), type: z.string(), /*value: z.any()*/ }).optional(),
      error: z.any().or(z.undefined()),
    },
  }, async (props) => {
    try {
      let fpl = props.workflow_id ? await fpprg.getFPL(props.workflow_id) : undefined
      if (props.workflow_id && fpl === undefined) throw new Error(`{"workflow_id":>"${props.workflow_id}" not found<}`)
      // assemble lookup with current steps in the graph
      const stepNodes = Object.fromEntries((fpl?.resolve() ?? []).map(step => [step.process.id, krg.getProcessNode(step.process.type)] as const))
      // Verify step is valid
      //  type should exist
      const processMetaNode = krg.getProcessNode(props.step.type)
      if (!processMetaNode) throw new Error(`{"step":{"type": >"${props.step.type}" is not from options< }}`)
      //  all specified derived_froms should be valid
      const errors: string[] = [] // we will collect any errors as we go to hopefully speed up the LLM getting it right
      for (const argName in props.step.derived_from) {
        const expectedArgType = processMetaNode.inputs[argName]
        if (expectedArgType === undefined) {
          errors.push(`{"step":{"derived_from":{">${argName}< is not expected": ...}, ...}}`)
          continue
        }
        const arg = props.step.derived_from[argName]
        const argNode = stepNodes[arg.step_id]
        if (!argNode) {
          errors.push(`{"step":{"derived_from":{"${argName}": {"step_id": >"${arg.step_id}" is not a valid step_id<, "type": "${expectedArgType.spec}" }, ...}, ...}}`)
          continue
        }
        if (argNode.output.spec !== expectedArgType.spec) {
          errors.push(`{"step":{"derived_from":{"${argName}":{"step_id": "${arg.step_id}", "type": >must be "${expectedArgType.spec}" got "${argNode.output.spec}"< }, ...}, ...}}`)
          continue
        }
      }
      //  all derived_froms should be specified
      for (const argName in processMetaNode.inputs) {
        if (props.step.derived_from[argName] === undefined) {
          errors.push(`{"step":{"derived_from": {"${argName}": {"step_id": >expected step_id<, "type": "${processMetaNode.inputs[argName].spec}" }, ...}, ...}}`)
          continue
        }
      }
      //  any values should be valid if a codec is defined
      if ('codec' in processMetaNode && props.step.value) {
        try {
          processMetaNode.codec.decode(props.step.value)
        } catch (e: any) {
          errors.push(e.toString())
        }
      } else {
        // ignore the value if provided when there is no codec
        props.step.value = undefined
      }
      if (errors.length > 0) throw new Error(errors.join('\n'))
      // actually register the step
      const process = await fpprg.resolveProcess({
        type: props.step.type,
        inputs: Object.fromEntries(Object.entries(props.step.derived_from).map(([key, value]) => [key, { id: value.step_id }])),
        data: props.step.value ? { type: processMetaNode.output.spec, value: props.step.value } : undefined,
      })
      if (!fpl) {
        fpl = await fpprg.upsertFPL(new FPL(process))
      } else {
        fpl = await fpprg.upsertFPL(fpl.extend(process))
      }
      const type = processMetaNode.output.spec
      // const output = await fpl.process.output()
      // if (!output) throw new Error('Failed to resolve output')
      // const { type, value } =  output.toJSON()
      const result = { workflow_id: fpl.id, step_id: fpl.process.id, type, /*value*/ }
      return {
        content: [{ type: 'text', text: JSON.stringify({ result }) }],
        structuredContent: { result },
      }
    } catch (exc) {
      const error = (exc as Error).message
      console.error(JSON.stringify({ expand: { input: props, error } }))
      return {
        content: [{ type: 'text', text: JSON.stringify({ error }) }],
        structuredContent: { error },
      }
    }
  })
  return server
})

export default handler(async (req, res) => {
  return await new Promise(async (resolve, reject) => {
    try {
      const transport = new StreamableHTTPServerTransport({
        sessionIdGenerator: undefined,
        enableJsonResponse: true,
      })
      res.on('close', () => {
        transport.close()
        resolve()
      })
      await server.connect(transport)
      await transport.handleRequest(req, res, req.body)
    } catch (e) {
      reject(e)
    }
  })
})
