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

const tools = {
  suggest_component: {
    description: 'Provide a list of suggested next steps to the user',
    parameters: {
      type: 'array',
      items: {
        type: 'string',
        description: 'A suggested component name from the available components',
      },
    },
    call: async (args: string) => {
      return ''
    },
  },
  extend_workflow: {
    description: 'Extend the workflow, returns a list of next steps which can now be added at this step.',
    parameters: {
      type: 'object',
      properties: {
        name: {
          type: 'string',
          description: 'A component name from the available components which will be appended to the current workflow',
        },
        data: {
          type: 'string',
          description: 'Instantiate the prompt with this data',
        },
      },
      required: ['name'],
    },
    call: async (args: string) => {
      const { name, data } = JSON.parse(args)
      return ''
    },
  },
}

const openai = cache('openai', () => new OpenAI({ apiKey: process.env.OPENAI_API_KEY }))
const assistant = cache('openai-assistant', async () => {
  if (!process.env.OPENAI_ASSISTANT_ID) {
    const assistant = await openai.beta.assistants.create({
      name: `playbook-gpt`,
      instructions: "You are a bot that helps the user construct workflows, you do not answer questions directly, instead you use the provided tools to construct workflows that answer the user's question.",
      model: "gpt-4-turbo-preview",
      tools: Object.entries(tools).map(([name, { call: _, ...tool }]) => ({ type: 'function', function: { name, ...tool } })),
    })
    console.log(assistant.id)
    return assistant
  } else {
    return await openai.beta.assistants.retrieve(process.env.OPENAI_ASSISTANT_ID)
  }
})

export const GPTAssistantThread3Create = API.post('/api/v1/chat/threads3')
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
        role: 'system',
        content: `You are a bot that helps the user construct workflows, you do not answer questions directly, instead you use the provided tools to construct workflows that answer the user's question.`,
      },
    })
    const thread = await openai.beta.threads.create()
    return msg.id
  })
  .build()

export const GPTAssistantThread3Message = API.post('/api/v1/chat/threads3/[thread_id]/messages')
  .query(z.object({ thread_id: z.string() }))
  .body(z.object({
    content: z.string(),
  }))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user || !process.env.OPENAI_API_KEY) throw new UnauthorizedError()
    if (req.method === 'HEAD') {
      res.status(200).end()
      return
    }
    const userMessage = await openai.beta.threads.messages.create(inputs.query.thread_id, {
      role: 'user',
      content: inputs.body.content,
    })
    const currentMessages = await openai.beta.threads.messages.list(inputs.query.thread_id, { before: userMessage.id, order: 'asc' })
    
    const assistant = await openai.beta.assistants.create({
      name: `playbook-gpt`,
      instructions: "You are a bot that helps the user construct workflows, you do not answer questions directly, instead you use the provided tools to construct workflows that answer the user's question.",
      model: "gpt-4-turbo-preview",
      tools: Object.entries(tools).map(([name, { call: _, ...tool }]) => ({ type: 'function', function: { name, ...tool } })),
    })

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

export const GPTAssistantThread3MessagesList = API.get('/api/v1/chat/threads3/[thread_id]/messages')
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

export const GPTAssistantThread3Delete = API.post('/api/v1/chat/threads3/[thread_id]/delete')
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
