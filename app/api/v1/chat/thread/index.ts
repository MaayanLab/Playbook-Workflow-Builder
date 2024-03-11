import { API } from '@/spec/api'
import { z } from 'zod'
import krg from '@/app/krg'
import * as dict from '@/utils/dict'
import db from '@/app/db'
import OpenAI from 'openai'
import type { ProcessMetaNode } from '@/spec/metanode'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { UnauthorizedError } from '@/spec/error'
import cache from '@/utils/global_cache'
import dedent from 'ts-dedent'

const openai = cache('openai', () => new OpenAI({ apiKey: process.env.OPENAI_API_KEY }))
const assistant = cache('openai-assistant', () => undefined /*openai.beta.assistants.create({
  name: `playbook-gpt`,
  instructions: dedent`
    You are a bot that helps the user construct workflows. You do not answer questions directly, instead you help choose from the applicable workflow steps in an effort to construct a workflow that answer the user's question.
    Queries will be of the form (note that steps in the current workflow are linked to the next step by id):
    {
      "message": "I'd like to get the expression of ACE2.",
      "choices": [{ "id": 1, "name": "Input[Gene]", ...]
    }
    or
    {
      "step": { "id": 1, "name": "Input[Gene]", "data": "ACE2" },
      "choices": [{ "id": 2, "name": "GetExpressionOfGene", "inputs": { "gene": 1 }, ... }, ...]
    }
    Your response should be of the form (choosing from the choices, suggestions should be ordered by most likely):
    {
      "agent_message": "You can get expression of ACE2 by ...",
      "suggestions": [{ "id": 2 }, ...]
    }
`,
  model: "gpt-4-turbo-preview",
})*/)

export const GPTAssistantThreadCreate = API.post('/api/v1/chat/threads')
  .query(z.object({}))
  .body(z.object({}))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user || !process.env.OPENAI_API_KEY) throw new UnauthorizedError()
    if (req.method === 'HEAD') {
      res.status(200).end()
      return
    }
    const msg = await db.objects.chat_message.create({
      data: {
        user: session.user.id,
        role: 'system',
        content: `You are a bot that helps the user construct workflows, you do not answer questions directly, instead you use the provided tools to construct workflows that answer the user's question.`,
      },
    })
    return msg.id
  })
  .build()

export const GPTAssistantThreadMessage = API.post('/api/v1/chat/threads/[thread_id]/messages')
  .query(z.object({ thread_id: z.string() }))
  .body(z.object({
    role: z.enum(['user']),
    content: z.string(),
  }))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user || !process.env.OPENAI_API_KEY) throw new UnauthorizedError()
    if (req.method === 'HEAD') {
      res.status(200).end()
      return
    }
    // get all messages, + the new one the user just sent
    const messages: {
      id: string;
      user: string;
      previous: string;
      role: string;
      content: string;
      created: Date;
    }[] = []
    let where: { id?: string, previous?: string } = { id: inputs.query.thread_id }
    while (true) {
      const message = await db.objects.chat_message.findUnique({ where })
      if (message === null) break
      else messages.push(message)
      where = { previous: message.id }
    }
    const previousMessage = messages[messages.length-1]
    messages.push(
      await db.objects.chat_message.create({ data: { user: session.user.id, role: inputs.body.role, content: inputs.body.content, previous: previousMessage.id } })
    )
    const userMessagePos = messages.length
    const workflow = messages.filter(msg => msg.role === 'component').map(msg => JSON.parse(msg.content))
    const currentProc = workflow.length > 0 ? workflow[workflow.length-1].name : ''
    const nextProcOptions = krg.getNextProcess(currentProc)
      .map(proc => ({ name: proc.spec, type: 'prompt' in proc ? 'prompt' : 'resolver', label: proc.meta.label, description: proc.meta.description, story: proc.story({}) }))

    let run = await openai.beta.threads.runs.create(inputs.query.thread_id, { assistant_id: (await assistant).id })
    while (true) {
      await new Promise<void>((resolve, reject) => {setTimeout(() => {resolve()}, 500)})
      run = await openai.beta.threads.runs.retrieve(inputs.query.thread_id, run.id)
      if (run.status === 'completed') {
        // send all new messages since the user's message to the user
        const newMessages = await openai.beta.threads.messages.list(inputs.query.thread_id, { after: userMessage.id, order: 'asc' })
        return { data: newMessages.data }
      } else if (run.status === 'requires_action' && run.required_action?.type === 'submit_tool_outputs') {
        // add tool results to the thread
        await openai.beta.threads.runs.submitToolOutputs(
          inputs.query.thread_id,
          run.id,
          {
            tool_outputs: await Promise.all(
              // run all the tools
              run.required_action.submit_tool_outputs.tool_calls.map(async (tool) => {
                const fn = tools[tool.function.name as keyof typeof tools].call
                const output = await fn(tool.function.arguments)
                return { tool_call_id: tool.id, output }
              })
            ),
          }
        )
      } else if (run.status === 'in_progress' || run.status === 'queued') {
        continue
      } else {
        return { error: 'An OpenAI error occurred, try again later' }
      }
    }
  })
  .build()

export const GPTAssistantThreadMessagesList = API.get('/api/v1/chat/threads/[thread_id]/messages')
  .query(z.object({ thread_id: z.string() }))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user || !process.env.OPENAI_API_KEY) throw new UnauthorizedError()
    if (req.method === 'HEAD') {
      res.status(200).end()
      return
    }
    const messages = await openai.beta.threads.messages.list(inputs.query.thread_id, { order: 'asc' })
    return messages.data.map(datum => ({
      role: datum.role,
      content: datum.content,
    }))
  })
  .build()

export const GPTAssistantThreadDelete = API.post('/api/v1/chat/threads/[thread_id]/delete')
  .query(z.object({ thread_id: z.string() }))
  .body(z.object({}))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user || !process.env.OPENAI_API_KEY) throw new UnauthorizedError()
    if (req.method === 'HEAD') {
      res.status(200).end()
      return
    }
    await openai.beta.threads.del(inputs.query.thread_id)
  })
  .build()
