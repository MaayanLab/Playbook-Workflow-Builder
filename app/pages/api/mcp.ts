import handler, { RouteHandler } from '@/utils/next-rest'
import { McpServer, ResourceTemplate } from '@modelcontextprotocol/sdk/server/mcp.js';
import { StreamableHTTPServerTransport } from '@modelcontextprotocol/sdk/server/streamableHttp.js';
import cache from '@/utils/global_cache';
import krg from '@/app/krg';
import * as dict from '@/utils/dict'
import { z } from 'zod'

const server = cache('mcp', () => {
  const server = new McpServer({
    name: 'Playbook Workflow Builder MCP Server',
    version: '1.0.0',
  })
  krg.getProcessNodes().forEach((metanode) => {
    server.registerTool(
      metanode.spec,
      {
        title: metanode.meta.label,
        description: metanode.meta.description,
        inputSchema: {
          ...('codec' in metanode ? {data: metanode.codec.zod ?? z.any()} : {}),
          inputs: z.object(dict.init(
            dict.items(metanode.inputs)
              .map(({ key, value }) => ({
                key,
                value: Array.isArray(value) ? (value[0].codec.zod ?? z.any()).array() : value.codec.zod ?? z.any(),
              }))
          )),
        },
        outputSchema: {
          result: (metanode.output.codec.zod ?? z.any()).or(z.undefined()),
          error: z.any().or(z.undefined()),
        },
      },
      async (props: any) => {
        console.warn(`MCP ${metanode.spec}: ${JSON.stringify(props)}`)
        try {
          const output = { result: await metanode.resolve({ ...props, notify(status) {} }) }
          return {
            content: [{ type: 'text', text: JSON.stringify(output) }],
            structuredContent: output as any
          }
        } catch (e) {
          console.error(e)
          const output = { error: (e as Error).toString() }
          return {
            content: [{ type: 'text', text: JSON.stringify(output) }],
            structuredContent: output
          }
        }
      }
    )
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
