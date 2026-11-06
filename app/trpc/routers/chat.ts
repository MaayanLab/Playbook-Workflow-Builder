import { z } from 'zod'
import OpenAI from 'openai'
import dedent from 'ts-dedent'
import { TRPCError } from '@trpc/server'
import db from '@/app/db'
import krg from '@/app/krg'
import fpprg from '@/app/fpprg'
import { uuid, type ThreadMessage } from '@/db'
import cache from '@/utils/global_cache'
import { authedProcedure, router } from '../init'

const openai = cache('openai', async () => {
  if (!process.env.OPENAI_API_KEY) {
    console.warn(`OPENAI_API_KEY not defined`)
    throw new TRPCError({ code: 'METHOD_NOT_SUPPORTED', message: 'OPENAI_API_KEY not defined' })
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

async function getThread(thread_id: string) {
  const thread = await db.selectFrom('thread').selectAll().where('id', '=', thread_id).executeTakeFirst()
  if (!thread) throw new TRPCError({ code: 'NOT_FOUND', message: thread_id })
  return thread
}

async function createMessage(message: Omit<ThreadMessage, 'id' | 'created' | 'feedback' | 'fpl'> & { fpl?: string | null }) {
  return await db.insertInto('thread_message')
    .values({ id: uuid(), ...message })
    .returningAll()
    .executeTakeFirstOrThrow()
}

async function updateMessage(id: string, data: Partial<Pick<ThreadMessage, 'fpl' | 'content'>>) {
  return await db.updateTable('thread_message')
    .set(data)
    .where('id', '=', id)
    .returningAll()
    .executeTakeFirstOrThrow()
}

// TODO: derive this from the mcp server definition somehow
const McpCallC = z.union([
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
])

export const chatRouter = router({
  available: authedProcedure
    .query(() => Boolean(process.env.OPENAI_API_KEY)),

  augment: authedProcedure
    .input(z.object({ description: z.string() }))
    .mutation(async ({ input }) => {
      if (!process.env.OPENAI_API_KEY) throw new TRPCError({ code: 'METHOD_NOT_SUPPORTED', message: 'OpenAI is not configured' })
      const response = await fetch('https://api.openai.com/v1/chat/completions', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json', Authorization: `Bearer ${process.env.OPENAI_API_KEY}` },
        body: JSON.stringify({
          model: 'gpt-4o',
          messages: [{ role: 'user', content: `Restructure this text so it is easier to read, do not change or omit any information:\n${input.description}` }],
          temperature: 0.7,
        }),
      })
      if (!response.ok) throw new TRPCError({ code: 'BAD_GATEWAY', message: await response.text() })
      const result = z.object({ choices: z.array(z.object({ message: z.object({ content: z.string() }) })) }).parse(await response.json())
      return result.choices[0].message.content.trim()
    }),

  create: authedProcedure
    .mutation(async ({ ctx }) => {
      const thread = await db.insertInto('thread')
        .values({ id: uuid(), user: ctx.session.user.id })
        .returning('id')
        .executeTakeFirstOrThrow()
      return thread.id
    }),

  messages: authedProcedure
    .input(z.object({ thread_id: z.string() }))
    .query(async ({ input }) => {
      const thread = await getThread(input.thread_id)
      const messages = await db.selectFrom('thread_message')
        .selectAll()
        .where('thread', '=', thread.id)
        .orderBy('created', 'asc')
        .execute()
      return { messages, fpl: messages[messages.length-1]?.fpl ?? null }
    }),

  send: authedProcedure
    .input(z.object({
      thread_id: z.string(),
      graph_id: z.string().optional(),
      node_id: z.string().optional(),
      message: z.string(),
    }))
    .subscription(async function* ({ input: { thread_id, message, ...input }, signal }) {
      const graph_id = input.graph_id !== 'start' ? input.graph_id : undefined
      const node_id = input.node_id !== 'start' ? input.node_id : undefined
      const thread = await getThread(thread_id)
      const currentMessages = await db.selectFrom('thread_message')
        .selectAll()
        .where('thread', '=', thread.id)
        .orderBy('created', 'asc')
        .execute()
      const newMessages: ThreadMessage[] = []
      const steps = currentMessages.flatMap(msg => msg.fpl ? [msg.fpl] : [])
      if (graph_id && !steps.includes(graph_id)) {
        const fpl = await fpprg.getFPL(graph_id)
        if (!fpl) throw new TRPCError({ code: 'NOT_FOUND', message: graph_id })
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
          newMessages.push(await createMessage({ thread: thread.id, ...newMessage }))
        }
      }
      const focusMessages: { role: string, content: string }[] = []
      if (node_id && graph_id !== node_id) {
        const fpl = await fpprg.getFPL(node_id)
        if (!fpl) throw new TRPCError({ code: 'NOT_FOUND', message: node_id })
        focusMessages.push(
          { role: 'developer', content: `User is currently looking at step ${fpl.process.id}` }
        )
      }
      newMessages.push(await createMessage({
        thread: thread.id,
        role: 'user',
        content: message,
        fpl: graph_id,
      }))
      yield { messages: newMessages, fpl: graph_id }
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
      }, { signal })
      let fpl: string | undefined
      const partialMessages: Record<number, ThreadMessage> = {}
      for await (const event of stream) {
        console.debug(JSON.stringify(event))
        if (event.type === 'response.output_item.added') {
          if (event.item.type === 'mcp_call') {
            const output = event.item
            partialMessages[event.output_index] = await createMessage({
              thread: thread.id,
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
            })
            newMessages.push(partialMessages[event.output_index])
            yield { messages: newMessages, fpl }
          } else if (event.item.type === 'message') {
            partialMessages[event.output_index] = await createMessage({
              thread: thread.id,
              role: 'assistant',
              content: event.item.content.flatMap(c => c.type === 'output_text' ? [c.text] : c.type === 'refusal' ? [c.refusal] : []).join(''),
              fpl,
            })
            newMessages.push(partialMessages[event.output_index])
            yield { messages: newMessages, fpl }
          }
        } else if (event.type === 'response.output_text.delta') {
          partialMessages[event.output_index] = await updateMessage(partialMessages[event.output_index].id, {
            fpl,
            content: partialMessages[event.output_index].content + event.delta,
          })
          yield { messages: newMessages, fpl }
        } else if (event.type === 'response.output_item.done') {
          if (event.item.type === 'mcp_call') {
            const output = event.item
            const outputParsed = McpCallC.safeParse({
              name: output.name,
              arguments: JSONTryParse(output.arguments),
              output: output.output ? JSONTryParse(output.output) : undefined,
            })
            fpl = outputParsed.data?.name === 'expand' ? outputParsed.data?.output.result.workflow_id
              : outputParsed.data?.name === 'options' ? outputParsed.data?.arguments.workflow_id
              : outputParsed.data?.name === 'view' ? outputParsed.data?.arguments.workflow_id
              : undefined
            if (event.output_index in partialMessages) {
              partialMessages[event.output_index] = await updateMessage(partialMessages[event.output_index].id, {
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
              })
              yield { messages: newMessages, fpl }
            }
          } else if (event.item.type === 'message') {
            partialMessages[event.output_index] = await updateMessage(partialMessages[event.output_index].id, {
              fpl,
              content: event.item.content.flatMap(c => c.type === 'output_text' ? [c.text] : c.type === 'refusal' ? [c.refusal] : []).join(''),
            })
            yield { messages: newMessages, fpl }
          }
        }
      }
    }),

  feedback: authedProcedure
    .input(z.object({ thread_id: z.string(), message_id: z.string(), feedback: z.string() }))
    .mutation(async ({ input }) => {
      await db.updateTable('thread_message')
        .set({ feedback: input.feedback })
        .where('id', '=', input.message_id)
        .where('thread', '=', input.thread_id)
        .execute()
      return null
    }),

  delete: authedProcedure
    .input(z.object({ thread_id: z.string() }))
    .mutation(async ({ input }) => {
      const thread = await getThread(input.thread_id)
      await db.deleteFrom('thread').where('id', '=', thread.id).execute()
      return null
    }),
})
