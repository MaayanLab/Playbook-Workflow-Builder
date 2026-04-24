import { API } from '@/spec/api'
import { z } from 'zod'
import db from '@/app/db'
import krg from '@/app/krg'
import OpenAI from 'openai'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { NotFoundError, UnauthorizedError, UnsupportedMethodError } from '@/spec/error'
import cache from '@/utils/global_cache'
import fpprg from '@/app/fpprg'

const openai = cache('openai', async () => {
  if (!process.env.OPENAI_API_KEY) {
    console.warn(`OPENAI_API_KEY not defined`)
    throw new UnsupportedMethodError()
  }
  return new OpenAI({ apiKey: process.env.OPENAI_API_KEY })
})

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
    const newMessages: { role: string, content: string }[] = []
    const steps = currentMessages.flatMap(msg => msg.fpl ? [msg.fpl] : [])
    if (inputs.body.graph_id && steps.includes(inputs.body.graph_id)) {
      const fpl = await fpprg.getFPL(inputs.body.graph_id)
      if (!fpl) throw new NotFoundError(inputs.body.graph_id)
      const graph_steps = fpl.resolve()
      for (const step of graph_steps.toReversed()) {
        if (steps.includes(step.id)) break
        const metaProcess = krg.getProcessNode(fpl.process.id)
        const newMessage = {
          thread: pwb_thread.id,
          fpl: step.id,
          role: 'developer',
          content: `triggered\n${JSON.stringify({
            function_call: {
              name: 'extend',
              arguments: {
                step: {
                  type: step.process.type,
                  derived_from: Object.fromEntries(Object.entries(step.process.inputs).map(([key, value]) => [key, { id: value.id }])),
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
        }
        newMessages.push(await db.objects.thread_message.create({ data: newMessage }))
      }
    }
    if (inputs.body.node_id && inputs.body.graph_id !== inputs.body.node_id) {
      const fpl = await fpprg.getFPL(inputs.body.node_id)
      if (!fpl) throw new NotFoundError(inputs.body.graph_id)
      newMessages.push(
        { role: 'developer', content: `User is currently looking at step ${fpl.process.id}` }
      )
    }
    const userMessage = await db.objects.thread_message.create({
      data: {
        thread: pwb_thread.id,
        role: 'user',
        content: inputs.body.message,
        fpl: inputs.body.graph_id,
      }
    })
    const response = await (await openai).responses.create({
      model: 'gpt-5-mini',
      tools: [
        {
          type: 'mcp',
          server_label: 'PWBMCP',
          server_url: `${process.env.PUBLIC_URL}/mcp`,
          require_approval: 'never',
        }
      ],
      instructions: 'You are an assistant that builds workflows with the playbook worfklow builder on behalf of the user using the PWBMCP server tools. Do NOT propose workflows without calling the PWB options tool. Do NOT interpret output. Only provide the constructed report for the user. The user can see when you `expand`, so do NOT include specific workflow/step IDs in the output.',
      input: [
        ...currentMessages,
        ...newMessages,
        userMessage,
      ].map(({ role, content }) => ({ role, content })) as OpenAI.Responses.ResponseInput,
    })
    let fpl: string | undefined
    for (const output of response.output) {
      if (output.type === 'mcp_call' && output.name === 'expand' && output.output) {
        try {
          const outputParsed = z.object({
            result: z.object({
              workflow_id: z.string(),
              step_id: z.string(),
              type: z.string(),
            }),
          }).parse(JSON.parse(output.output))
          fpl = outputParsed.result.workflow_id
          await db.objects.thread_message.create({ data: {
            thread: pwb_thread.id,
            fpl: outputParsed.result.workflow_id,
            role: 'developer',
            content: `${JSON.stringify({
              function_call: {
                name: 'extend',
                arguments: output.arguments,
              },
              function_call_output: {
                output: output.output,
              },
            })}`
          } })
        } catch (e: any) {}
      }
    }
    const assistantMessage = await db.objects.thread_message.create({
      data: {
        thread: pwb_thread.id,
        role: 'assistant',
        content: response.output_text,
        fpl,
      }
    })
    return {
      messages: [
        userMessage,
        assistantMessage,
      ],
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
    return { messages: pwb_thread_messages.filter(msg => msg.role !== 'developer'), fpl: pwb_thread_messages[pwb_thread_messages.length-1]?.fpl ?? null }
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
