import { API } from '@/spec/api'
import { z } from 'zod'
import db from '@/app/db'
import krg from '@/app/krg'
import OpenAI from 'openai'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { NotFoundError, UnauthorizedError, UnsupportedMethodError } from '@/spec/error'
import cache from '@/utils/global_cache'
import fpprg from '@/app/fpprg'
import dedent from 'ts-dedent'

const openai = cache('openai', async () => {
  if (!process.env.OPENAI_API_KEY) {
    console.warn(`OPENAI_API_KEY not defined`)
    throw new UnsupportedMethodError()
  }
  return new OpenAI({ apiKey: process.env.OPENAI_API_KEY })
})

function JSONTryParse(str: string) {
  try {
    return JSON.parse(str)
  } catch (e) {
    return undefined
  }
}

export const GPTAssistantCreate = API.post('/api/v1/chat')
  .query(z.object({}))
  .body(z.object({}))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    if (req.method === 'HEAD') {
      res.status(200).end()
      return
    }
    const pwb_thread = await db.objects.thread.create({ data: { user: session.user.id } })
    return pwb_thread.id
  })
  .build()

export const GPTAssistantMessage = API.post('/api/v1/chat/[thread_id]/messages')
  .query(z.object({
    thread_id: z.string(),
  }))
  .body(z.object({
    graph_id: z.string().optional(),
    node_id: z.string().optional(),
    message: z.string(),
  }))
  .call(async (inputs, req, res) => {
    if (inputs.body.graph_id === 'start') delete inputs.body.graph_id
    if (inputs.body.node_id === 'start') delete inputs.body.node_id
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    if (req.method === 'HEAD') {
      res.status(200).end()
      return
    }
    const pwb_thread = await db.objects.thread.findUnique({ where: { id: inputs.query.thread_id } })
    if (!pwb_thread) throw new NotFoundError(inputs.query.thread_id)
    const currentMessages = await db.objects.thread_message.findMany({
      where: { thread: pwb_thread.id },
      orderBy: { created: 'asc' },
    })
    const newMessages: typeof currentMessages = []
    const steps = currentMessages.flatMap(msg => msg.fpl ? [msg.fpl] : [])
    if (inputs.body.graph_id && !steps.includes(inputs.body.graph_id)) {
      const fpl = await fpprg.getFPL(inputs.body.graph_id)
      if (!fpl) throw new NotFoundError(inputs.body.graph_id)
      const graph_steps = fpl.resolve()
      const newMessageQueue: { fpl: string, role: string, content: string }[] = []
      for (const step of graph_steps.toReversed()) {
        if (steps.includes(step.id)) break
        const metaProcess = krg.getProcessNode(fpl.process.type)
        newMessageQueue.push({
          fpl: step.id,
          role: 'developer',
          content: `${JSON.stringify({
            function_call: {
              name: 'extend',
              arguments: {
                step: {
                  type: step.process.type,
                  arguments: Object.fromEntries(Object.entries(step.process.inputs).map(([key, value]) => [key, { id: value.id }])),
                  value: step.process.data,
                },
              },
            },
            function_call_output: {
              output: {
                result: {
                  workflow_id: step.id,
                  step_id: step.process.id,
                  output: { type: metaProcess.output.spec },
                },
              },
            },
          })}`
        })
      }
      for (const newMessage of newMessageQueue.toReversed()) {
        newMessages.push(await db.objects.thread_message.create({ data: { thread: pwb_thread.id, ...newMessage } }))
      }
    }
    const focusMessages: { role: string, content: string }[] = []
    if (inputs.body.node_id && inputs.body.graph_id !== inputs.body.node_id) {
      const fpl = await fpprg.getFPL(inputs.body.node_id)
      if (!fpl) throw new NotFoundError(inputs.body.node_id)
      focusMessages.push(
        { role: 'developer', content: `User is currently looking at step ${fpl.process.id}` }
      )
    }
    newMessages.push(await db.objects.thread_message.create({
      data: {
        thread: pwb_thread.id,
        role: 'user',
        content: inputs.body.message,
        fpl: inputs.body.graph_id,
      }
    }))
    const stream = await (await openai).responses.create({
      model: 'gpt-5-mini',
      tools: [
        {
          type: 'mcp',
          server_label: 'PWBMCP',
          server_url: `${process.env.PUBLIC_URL}/mcp`,
          require_approval: 'never',
        }
      ],
      instructions: dedent`
        You are an assistant that builds workflows with the playbook workflow builder on behalf of the user using the PWBMCP server tools.
        Check for applicable published workflows before building a new one, if so 'view' them and respond to the user.
        Otherwise, iteratively build the workflow with the user. See the MCP server instructions for more information.
        The user sees the workflow, you do not, so there is nothing for you to directly interpret, analyze, or provide for the user.
        The user can see when you 'view' or 'extend', so do NOT include specific workflow/step IDs in the output.
        Do NOT offer to provide results/exports, the user can do this from the UI.
        Do NOT assume capabilities, only what can be found with options for a given step is possible.
      `,
      input: [
        ...currentMessages,
        ...newMessages,
        ...focusMessages,
      ].map(({ role, content }) => ({ role, content })) as OpenAI.Responses.ResponseInput,
      stream: true,
    })
    let fpl: string | undefined
    const partialMessages: Record<number, any> = {}
    for await (const event of stream) {
      console.debug(JSON.stringify(event))
      if (event.type === 'response.output_item.added') {
        if (event.item.type === 'mcp_call') {
          const output = event.item
          partialMessages[event.output_index] = await db.objects.thread_message.create({
            data: {
              thread: pwb_thread.id,
              role: 'developer',
              content: `${JSON.stringify({
                function_call: {
                  name: output.name,
                  arguments: output.arguments,
                },
                function_call_output: {
                  output: output.output,
                },
              })}`
            }
          })
          newMessages.push(partialMessages[event.output_index])
        } else if (event.item.type === 'message') {
          partialMessages[event.output_index] = await db.objects.thread_message.create({
            data: {
              thread: pwb_thread.id,
              role: 'assistant',
              content: event.item.content.flatMap(c => c.type === 'output_text' ? [c.text] : c.type === 'refusal' ? [c.refusal] : []).join(''),
              fpl,
            }
          })
          newMessages.push(partialMessages[event.output_index])
        }
      }
      else if (event.type === 'response.output_text.delta') {
        partialMessages[event.output_index] = await db.objects.thread_message.update({
          where: { id: partialMessages[event.output_index].id },
          data: {
            fpl,
            content: partialMessages[event.output_index].content + event.delta,
          }
        })
      } else if (event.type === 'response.output_item.done') {
        if (event.item.type === 'mcp_call') {
          const output = event.item
          // TODO: derive this from the mcp server definition somehow
          const outputParsed = z.union([
            z.object({
              name: z.literal('options'),
              arguments: z.object({
                workflow_id: z.string().nullish().transform(val => val === null ? undefined : val),
                step_id: z.string().nullish(),
              })
            }).strip(),
            z.object({
              name: z.literal('expand'),
              output: z.object({
                result: z.object({
                  workflow_id: z.string(),
                  step_id: z.string(),
                  type: z.string(),
                }),
              })
            }).strip(),
            z.object({
              name: z.literal('view'),
              arguments: z.object({
                workflow_id: z.string().nullish().transform(val => val === null ? undefined : val),
              })
            }).strip(),
          ] ).safeParse({
            name: output.name,
            arguments: JSONTryParse(output.arguments),
            output: output.output ? JSONTryParse(output.output) : undefined,
          })
          fpl = outputParsed.data?.name === 'expand' ? outputParsed.data?.output.result.workflow_id
            : outputParsed.data?.name === 'options' ? outputParsed.data?.arguments.workflow_id
            : outputParsed.data?.name === 'view' ? outputParsed.data?.arguments.workflow_id
            : undefined
          if (event.output_index in partialMessages) {
            partialMessages[event.output_index] = await db.objects.thread_message.update({
              where: { id: partialMessages[event.output_index].id },
              data: {
                fpl,
                content: `${JSON.stringify({
                  function_call: {
                    name: output.name,
                    arguments: output.arguments,
                  },
                  function_call_output: {
                    output: output.output,
                  },
                })}`
              }
            })
          }
        } else if (event.item.type === 'message') {
          partialMessages[event.output_index] = await db.objects.thread_message.update({
            where: { id: partialMessages[event.output_index].id },
            data: {
              fpl,
              content: event.item.content.flatMap(c => c.type === 'output_text' ? [c.text] : c.type === 'refusal' ? [c.refusal] : []).join(''),
            }
          })
        }
      }
    }
    return {
      messages: newMessages,
      fpl,
    }
  })
  .build()

export const GPTAssistantMessageFeedback = API.post('/api/v1/chat/[thread_id]/messages/[message_id]/feedback')
  .query(z.object({ thread_id: z.string(), message_id: z.string() }))
  .body(z.string())
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    await db.objects.thread_message.update({
      where: { id: inputs.query.message_id, thread: inputs.query.thread_id },
      data: { feedback: inputs.body },
    })
    return null
  })
  .build()

export const GPTAssistantMessagesList = API.get('/api/v1/chat/[thread_id]/messages')
  .query(z.object({ thread_id: z.string() }))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    if (req.method === 'HEAD') {
      res.status(200).end()
      return
    }
    const pwb_thread = await db.objects.thread.findUnique({ where: { id: inputs.query.thread_id } })
    if (!pwb_thread) throw new NotFoundError(inputs.query.thread_id)
    const pwb_thread_messages = await db.objects.thread_message.findMany({ where: { thread: inputs.query.thread_id }, orderBy: { created: 'asc' } })
    return { messages: pwb_thread_messages, fpl: pwb_thread_messages[pwb_thread_messages.length-1]?.fpl ?? null }
  })
  .build()

export const GPTAssistantDelete = API.post('/api/v1/chat/[thread_id]/delete')
  .query(z.object({ thread_id: z.string() }))
  .body(z.object({}))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    if (req.method === 'HEAD') {
      res.status(200).end()
      return
    }
    const pwb_thread = await db.objects.thread.findUnique({ where: { id: inputs.query.thread_id } })
    if (!pwb_thread) throw new NotFoundError(inputs.query.thread_id)
    await db.objects.thread.delete({ where: { id: inputs.query.thread_id } })
    return null
  })
  .build()
