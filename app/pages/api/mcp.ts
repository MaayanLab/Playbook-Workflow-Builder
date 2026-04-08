import handler, { RouteHandler } from '@/utils/next-rest'
import { McpServer, ResourceTemplate } from '@modelcontextprotocol/sdk/server/mcp.js';
import { StreamableHTTPServerTransport } from '@modelcontextprotocol/sdk/server/streamableHttp.js';
import cache from '@/utils/global_cache';
import krg from '@/app/krg';
import * as dict from '@/utils/dict'
import { z } from 'zod'
import dedent from 'ts-dedent';
import * as array from '@/utils/array'
import { ProcessMetaNode } from '@/spec/metanode';
import fpprg from '@/app/fpprg';
import { Metapath } from '@/app/fragments/metapath';
import { FPL } from '@/core/FPPRG';

function count(start: number) {
  const ctx = { id: start }
  return {
    next: () => ctx.id++
  }
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

      At each step, use the \`expand\` tool with the current workflow and choose from the next possible options with confirmation from the user.

      ## Important Notes

      PWB has a network of possibilities. The next possible option is not always what the user wants as an end product, but it might be a means to get there.
    `,
  })
  server.registerTool('options', {
    title: 'Options to start or expand a workflow with',
    inputSchema: {
      workflow_id: z.string().optional().describe('The complete workflow we are building'),
      step_id: z.string().optional().describe('A specific step in the workflow we wish to extend from'),
    },
    outputSchema: {
      result: z.object({ name: z.string(), inputs: z.record(z.string(), z.object({ id: z.string().describe('output of previous step_id') })) }).optional(),
      error: z.any().or(z.undefined()),
    },
  }, async (props) => {
    try {
      const fpl = props.workflow_id ? await fpprg.getFPL(props.workflow_id) : undefined
      if (props.workflow_id && !fpl) throw new Error('Workflow not found')
      const head_id = props.step_id ? props.step_id : props.workflow_id
      const workflow = fpl?.resolve() ?? []
      let head
      for (const item of workflow ?? []) {
        if (item.id === head_id) {
          head = item
          break
        }
      }
      if (props.step_id && !head) throw new Error('Step not found')
      
      const selections: Record<string, { process: typeof workflow[0], processNode: ProcessMetaNode }> = {}
      for (const item of workflow) {
        if (item === undefined) continue
        // add this to the selections
        selections[item.id] = { process: item, processNode: krg.getProcessNode(item.process.type) }
        // if a selection previously registered is a parent of this selection, remove it from selections
        dict.values(await item.process.inputs__outputs()).forEach(k => {
          if (!k) return
          if (k.id in selections) delete selections[k.id]
        })
      }
      const choices = krg.getNextProcess(head?.process.type)
        .filter(item => item.meta.hidden !== true)
        .filter(item => array.all(
          dict.values(item.inputs).map((value) => {
            if (Array.isArray(value)) {
              return dict.values(selections).filter(selection => selection.processNode.output.spec === value[0].spec).length > 1
            } else {
              return dict.values(selections).filter(selection => selection.processNode.output.spec === value.spec).length >= 1
            }
          })
        ))
        .map(item => {
          const inputs: Record<string, { id: string }> = {}
          dict.items(item.inputs).forEach(({ key: arg, value: input }) => {
            if (Array.isArray(input)) {
              dict.values(selections)
                .filter(selection => selection.processNode.output.spec === input[0].spec)
                .forEach((selection, i) => {
                  inputs[`${arg}:${i}`] = { id: selection.process.id }
                })
            } else {
              const head = { process: workflow[workflow.length-1] }
              const relevantSelections = dict.filter(selections, ({ value: selection }) => selection.processNode.output.spec === input.spec)
              const selection = head.process.id in relevantSelections ? head : array.ensureOne(dict.values(relevantSelections))
              inputs[arg] = { id: selection.process.id }
            }
          })
          return {
            name: item.spec, inputs,
            type: 'prompt' in item ? 'prompt' : 'resolver',
            label: item.meta.label,
            description: item.story({}) ?? item.meta.description,
          }
        })
      return {
        content: [{ type: 'text', text: JSON.stringify({ result: choices }) }],
        structuredContent: { result: choices },
      }
    } catch (error) {
      return {
        content: [{ type: 'text', text: JSON.stringify({ error }) }],
        structuredContent: { error },
      }
    }
  })
  server.registerTool('expand', {
    title: 'Expand Workflow',
    inputSchema: {
      workflow_id: z.string().or(z.literal('start')).describe('The complete workflow we are building off of'),
      step_id: z.string().optional().describe('A specific step in the workflow we wish to extend from'),
      step: z.object({
        name: z.string(),
        value: z.string().optional(),
        inputs: z.record(z.string(), z.object({ id: z.string() })),
      }).describe('What we do at this step'),
    },
    outputSchema: {
      result: z.object({ workflow_id: z.string(), step_id: z.string(), output: z.string() }).optional(),
      error: z.any().or(z.undefined()),
    },
  }, async (props) => {
    try {
      const process = await fpprg.resolveProcess({
        type: props.step.name,
        inputs: props.step.inputs,
      })
      let fpl: FPL
      if (props.workflow_id === 'start') {
        fpl = await fpprg.upsertFPL(new FPL(process))
      } else {
        const old_fpl = await fpprg.getFPL(props.workflow_id)
        if (old_fpl === undefined) throw new Error('Workflow not found')
        fpl = await fpprg.upsertFPL(old_fpl.extend(process))
      }
      const output = await fpl.process.output()
      if (!output) throw new Error('Failed to resolve output')
      const { type, value } =  output?.toJSON()
      return {
        content: [{ type: 'text', text: JSON.stringify({ workflow_id: fpl.id, step_id: fpl.process.id, output: { type, value } }) }],
        structuredContent: { workflow_id: fpl.id, step_id: fpl.process.id },
      }
    } catch (error) {
      return {
        content: [{ type: 'text', text: JSON.stringify({ error }) }],
        structuredContent: { error },
      }
    }
  })
  return server
})

export default handler(async (req, res) => {
  const transport = new StreamableHTTPServerTransport({
    sessionIdGenerator: undefined,
    enableJsonResponse: true,
  })
  res.on('close', () => {transport.close()})
  await server.connect(transport)
  await transport.handleRequest(req, res, req.body)
})
