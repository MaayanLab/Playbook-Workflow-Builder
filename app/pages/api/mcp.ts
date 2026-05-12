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
import { PublicPlaybooks, PublicUserPlaybooks } from '@/app/api/server'
import { fpl_expand } from '@/core/common';

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

      Before creating a workflow, check to see if published workflows can be found with 'search_published' which already does what the user wants, otherwise proceed with the other tools.
      When building a workflow, use the 'options' tool with the most recent workflow and relevant step and choose from the next possible options using the 'expand' tool.  
      To start a workflow, use the 'options' without specifying any 'workflow_id' or 'step_id', and subsequently 'expand' with the relevant step without a 'workflow_id' or 'step_id'.  
      Step 'type' specified in 'expand' MUST appear in 'options' from the same workflow/step.
      Always use the latest workflow id, re-use prior step_ids when applicable.
      The workflow_id can be used to build a link for the user with: https://playbook-workflow-builder.cloud/report/{workflow_id}

      ## Example

      \`\`\`
      {"type": "user_message", "content": "I'd like to get the expression of ACE2."}
      {"type": "function_call", "name": "search_published", "arguments": {""search":"Gene expression"}}
      {"type": "function_call_output", "output": {"result":[]}}
      {"type": "function_call", "name": "options", "arguments": {}}
      {"type": "function_call_output", "output": {"result":[{"type":"Input[Gene]","arguments":{},"label":"Gene Input","description":"The workflow starts with selecting a search term."}]}}
      {"type": "user_message", "content": "I'd like to get the expression of ACE2."}
      {"type": "function_call", "name": "expand", "arguments": {"step":{"type":"Input[Gene]","arguments":{},"value":"ACE2"}}}
      {"type": "function_call_output", "output": {"result":{"workflow_id": "38b60ad6", "step_id": "2b976512", "type": "Gene"}}
      {"type": "function_call", "name": "options", "arguments": {"workflow_id":"38b60ad6", "step_id": "2b976512"}}
      {"type": "function_call_output", "output": {"result":[{"type":"GetTissueExpressionOfGene","arguments":{"gene": {"type":"Gene"}},"label":"Get expression of the gene", ...}]}}
      {"type": "function_call", "name": "expand", "arguments": {"workflow_id":"38b60ad6", "step":{"type":"GetTissueExpressionOfGene","arguments":{"gene": {"step_id": "2b976512"}}}}}
      {"type": "function_call_output", "output": {"result":{"workflow_id": "811b99c8", "step_id": "965ac0ed", "type": "TissueGeneExpression"}}
      {"type": "agent_message", "content": "I've created the report for you at <https://playbook-workflow-builder.cloud/report/811b99c8>, let me know what else you'd like to do."}
      \`\`\`

      ## Important Notes

      PWB has a network of possibilities. The next possible option is not always what the user wants as an end product, but it might be a means to get there.
    `,
  })
  server.registerTool('options', {
    title: 'Options',
    description: 'Next step types compatible with the current workflow step. See instructions for more information and an example.',
    inputSchema: {
      workflow_id: z.string().nullish().describe('Most recent workflow id (or otherwise null if no workflow yet present)'),
      step_id: z.string().nullish().describe('Step id in the workflow we intend to expand from (or otherwise null if no workflow yet present)'),
    },
    outputSchema: {
      result: z.object({
        type: z.string().describe('type to be used in extend'),
        label: z.string(),
        description: z.string(),
        arguments: z.record(z.string(), z.object({ type: z.string(), description: z.string() }).describe('output types this step allows'))
      }).array().optional(),
      error: z.any().or(z.undefined()),
    },
  }, async (props) => {
    try {
      console.debug(JSON.stringify({ options: { input: props } }))
      if (props.workflow_id === 'start') props.workflow_id = undefined
      if (props.step_id === 'start') props.step_id = undefined
      const fpl = props.workflow_id ? await fpprg.getFPL(props.workflow_id) : undefined
      if (props.workflow_id && !fpl) throw new Error(`{"workflow_id":>"${props.workflow_id}" not found<}`)
      const workflow = fpl?.resolve() ?? []
      let heads: FPL[] = []
      if (props.step_id) {
        for (const item of workflow ?? []) {
          if (item.process.id === props.step_id) {
            heads = [item]
            break
          }
        }
        if (heads.length === 0) throw new Error(`{"step_id":>"${props.step_id}" not found<}`)
      } else if (workflow.length > 0) {
        heads = [workflow[workflow.length-1]]
      }
      const { items } = options({ krg, heads, workflow })
      const result = items.map(proc => ({
        type: proc.spec,
        arguments: dict.init(dict.items(proc.inputs).map(({ key, value }) => ({
          key,
          value: Array.isArray(value) ? {
            type: value[0].spec, description: value[0].meta.description
          } : {
            type: value.spec, description: value.meta.description
          }
        }))),
        label: proc.meta.label,
        description: proc.story({}).abstract ?? proc.meta.description,
      }))
      console.debug(JSON.stringify({ options: { output: result } }))
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
    title: 'Create or Expand Workflow',
    description: 'Create new workflow or build off of current workflow. Use only after calling options first. See instructions for more information and an example.',
    inputSchema: {
      workflow_id: z.string().nullish().describe('Current workflow id (null for new workflow)'),
      step: z.object({
        type: z.string().describe('step type as seen in options'),
        arguments: z.record(z.string(), z.object({ step_id: z.string().describe('step_id which has a type the same as in options for this argument') }).strip()).describe('all arguments from options must be specified from a prior step').default({}),
        value: z.string().nullish().describe('input value for user when arguments are empty'),
      }).strict(),
    },
    outputSchema: {
      result: z.object({ workflow_id: z.string(), step_id: z.string(), type: z.string(), description: z.string() /*value: z.any()*/ }).optional(),
      error: z.any().or(z.undefined()),
    },
  }, async (props) => {
    try {
      console.debug(JSON.stringify({ options: { expand: props } }))
      let fpl = props.workflow_id ? await fpprg.getFPL(props.workflow_id) : undefined
      if (props.workflow_id && fpl === undefined) throw new Error(`{"workflow_id":>"${props.workflow_id}" not found<}`)
      // assemble lookup with current steps in the graph
      const stepNodes = Object.fromEntries((fpl?.resolve() ?? []).map(step => [step.process.id, krg.getProcessNode(step.process.type)] as const))
      // Verify step is valid
      //  type should exist
      const processMetaNode = krg.getProcessNode(props.step.type)
      if (!processMetaNode) throw new Error(`{"step":{"type": >"${props.step.type}" is not from options< }}`)
      //  all specified argumentss should be valid
      const errors: string[] = [] // we will collect any errors as we go to hopefully speed up the LLM getting it right
      for (const argName in props.step.arguments) {
        const expectedArgType = processMetaNode.inputs[argName]
        if (expectedArgType === undefined) {
          errors.push(`{"step":{"arguments":{">${argName}< is not expected": ...}, ...}}`)
          continue
        }
        const arg = props.step.arguments[argName]
        const argNode = stepNodes[arg.step_id]
        if (!argNode) {
          errors.push(`{"step":{"arguments":{"${argName}": {"step_id": >"${arg.step_id}" is not a valid step_id<, "type": "${expectedArgType.spec}" }, ...}, ...}}`)
          continue
        }
        if (argNode.output.spec !== expectedArgType.spec) {
          errors.push(`{"step":{"arguments":{"${argName}":{"step_id": "${arg.step_id}", "type": >must be "${expectedArgType.spec}" got "${argNode.output.spec}"< }, ...}, ...}}`)
          continue
        }
      }
      //  all argumentss should be specified
      for (const argName in processMetaNode.inputs) {
        if (props.step.arguments[argName] === undefined) {
          errors.push(`{"step":{"arguments": {"${argName}": {"step_id": >expected step_id<, "type": "${processMetaNode.inputs[argName].spec}" }, ...}, ...}}`)
          continue
        }
      }
      //  any values should be valid if a codec is defined
      if ('codec' in processMetaNode && props.step.value) {
        try {
          processMetaNode.codec.decode(props.step.value)
        } catch (e: any) {
          errors.push(`{"step":{"value": >${e.toString()}<, ...}, ...}`)
        }
      } else {
        // ignore the value if provided when there is no codec
        props.step.value = undefined
      }
      if (errors.length > 0) throw new Error(errors.join('\n'))
      // actually register the step
      const process = await fpprg.resolveProcess({
        type: props.step.type,
        inputs: Object.fromEntries(Object.entries(props.step.arguments).map(([key, value]) => [key, { id: value.step_id }])),
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
      const description = processMetaNode.story({
        data: process.data,
        inputs: await process.inputs__outputs(),
      }).abstract ?? processMetaNode.meta.description ?? ''
      const result = { workflow_id: fpl.id, step_id: fpl.process.id, type, description, /*value*/ }
      console.debug(JSON.stringify({ expand: { output: result } }))
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
  server.registerTool('search_published', {
    title: 'Search Published Workflow Apps',
    description: 'Find a published workflow that already does what the user wants. Use view to load the workflow and get more details.',
    inputSchema: {
      search: z.string().describe('A search that will identify relevant workflows, structure of the workflow like tools, databases, input/output types are relevant, specific biological entities are NOT.'),
    },
    outputSchema: {
      result: z.object({ workflow_id: z.string(), title: z.string(), description: z.string().optional(), dataSources: z.string().optional(), inputs: z.string().optional(), outputs: z.string().optional() }).array().optional(),
      error: z.any().or(z.undefined()),
    },
  }, async (props) => {
    try {
      const matchingPlaybooks = [
        ...(await PublicPlaybooks.call({ query: { search: props.search, limit: 10 }, body: undefined }, undefined as any, undefined as any)), // none of these undefined params matter, maybe in the future we can fix the typing..
        ...(await PublicUserPlaybooks.call({ query: { search: props.search, limit: 10 }, body: undefined }, undefined as any, undefined as any)), // none of these undefined params matter, maybe in the future we can fix the typing..
      ]
      const result = matchingPlaybooks.map((p) => ({ workflow_id: 'playbook' in p ? p.playbook : p.id, title: 'label' in p ? p.label : 'title' in p ? p.title : 'Untitled', description: p.description, inputs: p.inputs, outputs: p.outputs }))
      console.debug(JSON.stringify({ search_published: { output: result } }))
      return {
        content: [{ type: 'text', text: JSON.stringify({ result }) }],
        structuredContent: { result },
      }
    } catch (exc) {
      const error = (exc as Error).message
      console.error(JSON.stringify({ search_published: { input: props, error } }))
      return {
        content: [{ type: 'text', text: JSON.stringify({ error }) }],
        structuredContent: { error },
      }
    }
  })
  server.registerTool('view', {
    title: 'View a workflow by id',
    description: 'Load and view a workflow by its id.',
    inputSchema: {
      workflow_id: z.string(),
    },
    outputSchema: {
      result: z.object({
        workflow_id: z.string(),
        step_id: z.string(),
        arguments: z.record(z.string(), z.object({ step_id: z.string() })),
        type: z.string(),
        description: z.string(),
      }).array().optional(),
      error: z.any().or(z.undefined()),
    },
  }, async (props) => {
    try {
      const fpl = await fpprg.getFPL(props.workflow_id)
      if (!fpl) throw new Error(`{"workflow_id":>"${props.workflow_id}" not found<}`)
      const { fullFPL, processLookup } = await fpl_expand({ krg, fpl })
      const result = fullFPL
        .map((step) => {
          const { metanode, story } = processLookup[step.process.id]
          return {
            workflow_id: step.id,
            step_id: step.process.id,
            arguments: dict.init(dict.items(step.process.inputs).map(({ key, value }) => ({ key, value: { step_id: value.id } }))),
            type: metanode.output.spec,
            description: story.abstract ?? metanode.meta.description ?? '',
          }
        })
      console.debug(JSON.stringify({ view: { output: result } }))
      return {
        content: [{ type: 'text', text: JSON.stringify({ result }) }],
        structuredContent: { result },
      }
    } catch (exc) {
      const error = (exc as Error).message
      console.error(JSON.stringify({ view: { input: props, error } }))
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
