import { z } from 'zod'
import OpenAI from 'openai'
import dedent from 'ts-dedent'
import { TRPCError } from '@trpc/server'
import db from '@/app/db'
import { getKRG } from '@/app/server-krg'
import fpprg from '@/app/fpprg'
import { isPrompt, processArguments, reportLink, stepDescription, summarize } from '@/app/mcp-util'
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

async function createMessage(message: Omit<ThreadMessage, 'id' | 'created' | 'feedback' | 'fpl' | 'response_item'> & { fpl?: string | null, response_item?: string | null }) {
  return await db.insertInto('thread_message')
    .values({ id: uuid(), ...message })
    .returningAll()
    .executeTakeFirstOrThrow()
}

async function updateMessage(id: string, data: Partial<Pick<ThreadMessage, 'fpl' | 'content' | 'response_item'>>) {
  return await db.updateTable('thread_message')
    .set(data)
    .where('id', '=', id)
    .returningAll()
    .executeTakeFirstOrThrow()
}

type ResponseInputItem = OpenAI.Responses.ResponseInputItem

/** The errors which tRPC's subscription link reconnects after (see `retryableRpcCodes`) */
const retriedCodes: TRPCError['code'][] = ['INTERNAL_SERVER_ERROR', 'BAD_GATEWAY', 'SERVICE_UNAVAILABLE', 'GATEWAY_TIMEOUT']

function inputMessage(role: 'user' | 'developer' | 'system', content: string): ResponseInputItem {
  return { type: 'message', role, content: [{ type: 'input_text', text: content }] }
}

/**
 * Prefer lossless typed items, falling back to the projection used by old threads.
 *  Except the model's own output items (which have an `id`, e.g. its tool calls): OpenAI refuses
 *  them without the reasoning items they were produced with, which aren't kept, so they're given
 *  as their projection (their content describes them).
 */
function messageToResponseInput(message: ThreadMessage): ResponseInputItem {
  if (message.response_item) {
    const item = JSONTryParse(message.response_item)
    if (item && typeof item === 'object' && typeof item.type === 'string' && !('id' in item)) return item as ResponseInputItem
  }
  const role = message.role === 'assistant' ? 'assistant'
    : message.role === 'user' ? 'user'
    : 'developer'
  return { role, content: message.content }
}

/** Keep provider-specific context on the server; clients only receive the UI projection. */
function publicMessage({ response_item: _responseItem, ...message }: ThreadMessage) {
  return message
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
    name: z.union([z.literal('plan'), z.literal('set_value')]),
    output: z.object({
      result: z.object({
        workflow_id: z.string(),
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
      return { messages: messages.map(publicMessage), fpl: messages[messages.length-1]?.fpl ?? null }
    }),

  send: authedProcedure
    .input(z.object({
      thread_id: z.string(),
      graph_id: z.string().optional(),
      node_id: z.string().optional(),
      message: z.string(),
    }))
    .subscription(async function* ({ input: { thread_id, message, ...input }, signal }) {
      /** Save any streamed text which hasn't been saved yet (when the stream ends early) */
      let flushPartials = async () => {}
      // a subscription which fails with e.g. INTERNAL_SERVER_ERROR is reconnected (so re-sent) by the
      //  client (tRPC's httpSubscriptionLink): sending isn't safe to repeat (the message would be
      //  added & answered again), so a failure is reported as an error which isn't retried
      try {
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
          const krg = await getKRG()
          const graph_steps = fpl.resolve()
          const newMessageQueue: { fpl: string, role: string, content: string }[] = []
          for (const step of graph_steps.toReversed()) {
            if (steps.includes(step.id)) break
            // a step added in the UI is recorded as the `expand` call (see app/mcp.ts) which would have
            //  added it, so the model sees the user's decisions as if they were its own tool calls
            const metaProcess = krg.getProcessNode(step.process.type)
            let value: unknown = step.process.data?.value
            if (value !== undefined) {
              try { if (metaProcess && isPrompt(metaProcess)) value = metaProcess.codec.decode(value as any) } catch {}
              value = summarize(value)
            }
            const args = {
              ...(step.parent ? { workflow_id: step.parent.id } : {}),
              step: {
                type: step.process.type,
                arguments: processArguments(step.process),
                ...(value !== undefined ? { value } : {}),
              },
            }
            const result = {
              workflow_id: step.id,
              parent_workflow_id: step.parent?.id ?? null,
              link: reportLink(step.id),
              step_id: step.process.id,
              step_type: step.process.type,
              type: metaProcess?.output.spec ?? '?',
              data_type: metaProcess?.output.spec ?? '?',
              label: metaProcess?.meta.label ?? step.process.type,
              description: await stepDescription({ krg, fpprg, process: step.process }),
            }
            newMessageQueue.push({
              fpl: step.id,
              role: 'developer',
              content: JSON.stringify({
                function_call: { name: 'expand', arguments: JSON.stringify(args) },
                function_call_output: { output: JSON.stringify({ result }) },
              }),
            })
          }
          // (no response_item: the model is given the content, exactly like the projection of its own tool calls)
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
          response_item: JSON.stringify(inputMessage('user', message)),
        }))
        yield { messages: newMessages.map(publicMessage), fpl: graph_id }
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
            Check for applicable published workflows before building a new one, if so 'view' them (and use 'set_value' to run them on the user's inputs) and respond to the user.
            Otherwise, iteratively build the workflow with the user. See the MCP server instructions for more information.
            The user sees the workflow, you do not, so there is nothing for you to directly interpret, analyze, or provide for the user.
            The user can see when you 'view', 'expand', 'plan' or 'set_value', so do NOT include specific workflow/step IDs in the output.
            Do NOT offer to provide results/exports, the user can do this from the UI.
            Do NOT assume capabilities, only what can be found with options for a given step is possible.
          `,
          input: [
            ...currentMessages.map(messageToResponseInput),
            ...newMessages.map(messageToResponseInput),
            ...focusMessages.map(({ role, content }) => inputMessage(role as 'developer', content)),
          ],
          stream: true,
        }, { signal })
        let fpl: string | undefined
        const partialMessages: Record<number, ThreadMessage> = {}
        /** Replace a message (in what's sent to the client too, which holds the objects as they were added) */
        const setPartial = (index: number, message: ThreadMessage) => {
          partialMessages[index] = message
          const i = newMessages.findIndex(({ id }) => id === message.id)
          if (i >= 0) newMessages[i] = message
          else newMessages.push(message)
        }
        /** A tool call as the client renders it (app/fragments/chat/message.tsx) & the model is given it */
        const mcpCallContent = (item: { name: string, arguments: string, output?: string | null, error?: unknown }) => JSON.stringify({
          function_call: { name: item.name, arguments: item.arguments },
          function_call_output: { output: item.output, ...(item.error ? { error: item.error } : {}) },
        })
        const textContent = (item: { content: ({ type: 'output_text', text: string } | { type: 'refusal', refusal: string })[] }) =>
          item.content.flatMap(c => c.type === 'output_text' ? [c.text] : c.type === 'refusal' ? [c.refusal] : []).join('')
        // text deltas are sent to the client as they arrive, but saved at most this often (& when done)
        const SAVE_INTERVAL_MS = 1000
        const lastSaved: Record<number, number> = {}
        /** messages (by output index) with streamed text that isn't saved yet */
        const unsaved = new Set<number>()
        flushPartials = async () => {
          for (const index of unsaved) {
            const partial = partialMessages[index]
            try { await updateMessage(partial.id, { fpl: partial.fpl, content: partial.content }) } catch {}
          }
          unsaved.clear()
        }
        for await (const event of stream) {
          console.debug(JSON.stringify(event))
          if (event.type === 'response.output_item.added') {
            if (event.item.type === 'mcp_call') {
              setPartial(event.output_index, await createMessage({
                thread: thread.id,
                role: 'developer',
                content: mcpCallContent(event.item),
              }))
              yield { messages: newMessages.map(publicMessage), fpl }
            } else if (event.item.type === 'message') {
              setPartial(event.output_index, await createMessage({
                thread: thread.id,
                role: 'assistant',
                content: textContent(event.item),
                fpl,
              }))
              lastSaved[event.output_index] = Date.now()
              yield { messages: newMessages.map(publicMessage), fpl }
            }
          } else if (event.type === 'response.output_text.delta') {
            const partial = partialMessages[event.output_index]
            if (!partial) continue
            const content = partial.content + event.delta
            if (Date.now() - (lastSaved[event.output_index] ?? 0) >= SAVE_INTERVAL_MS) {
              setPartial(event.output_index, await updateMessage(partial.id, { fpl, content }))
              lastSaved[event.output_index] = Date.now()
              unsaved.delete(event.output_index)
            } else {
              setPartial(event.output_index, { ...partial, fpl: fpl ?? null, content })
              unsaved.add(event.output_index)
            }
            yield { messages: newMessages.map(publicMessage), fpl }
          } else if (event.type === 'response.output_item.done') {
            if (event.item.type === 'mcp_call') {
              const output = event.item
              const outputParsed = McpCallC.safeParse({
                name: output.name,
                arguments: JSONTryParse(output.arguments),
                output: output.output ? JSONTryParse(output.output) : undefined,
              })
              // keep the current workflow for calls which don't change it (search, errors, dry runs, ...)
              const call = outputParsed.data
              fpl = call === undefined ? fpl
                : 'output' in call ? call.output.result.workflow_id
                : call.arguments.workflow_id ?? fpl
              if (event.output_index in partialMessages) {
                setPartial(event.output_index, await updateMessage(partialMessages[event.output_index].id, {
                  fpl,
                  response_item: JSON.stringify(event.item),
                  content: mcpCallContent(output),
                }))
                yield { messages: newMessages.map(publicMessage), fpl }
              }
            } else if (event.item.type === 'message') {
              if (!(event.output_index in partialMessages)) continue
              setPartial(event.output_index, await updateMessage(partialMessages[event.output_index].id, {
                fpl,
                response_item: JSON.stringify(event.item),
                content: textContent(event.item),
              }))
              unsaved.delete(event.output_index)
              yield { messages: newMessages.map(publicMessage), fpl }
            }
          }
        }
      } catch (error) {
        if (signal?.aborted) return
        if (error instanceof TRPCError && !retriedCodes.includes(error.code)) throw error
        console.error('[chat.send]', error)
        throw new TRPCError({ code: 'UNPROCESSABLE_CONTENT', message: error instanceof Error ? error.message : String(error), cause: error })
      } finally {
        // (also when the client disconnects, which ends the generator without an error)
        await flushPartials()
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
