/**
 * Procedures used by CAVATICA workers to serve requests (see @/app/extensions/cavatica/tunnel)
 *
 * Workers authenticate with `Authorization: Token <key>`, the CAVATICA API key of
 *  the user who launched the session (which the worker was started with).
 */
import { z } from 'zod'
import { timingSafeEqual } from 'crypto'
import { TRPCError } from '@trpc/server'
import db from '@/app/db'
import { connectWorker, respond, type WorkerEvent } from '@/app/extensions/cavatica/tunnel'
import { publicProcedure, router, type Context } from '../init'

function sameSecret(a: string, b: string) {
  const [bufA, bufB] = [Buffer.from(a), Buffer.from(b)]
  return bufA.length === bufB.length && timingSafeEqual(bufA, bufB)
}

async function authorizeWorker(ctx: Context, session_id: string) {
  const token = /^Token (.+)$/.exec(ctx.req?.headers.authorization ?? '')?.[1]
  const session = await db.selectFrom('proxy_session')
    .innerJoin('user_integrations', 'user_integrations.id', 'proxy_session.user')
    .select('user_integrations.cavatica_api_key')
    .where('proxy_session.id', '=', session_id)
    .executeTakeFirst()
  if (!token || !session?.cavatica_api_key || !sameSecret(token, session.cavatica_api_key)) {
    throw new TRPCError({ code: 'UNAUTHORIZED' })
  }
}

const SessionInput = z.object({ session_id: z.string() })

export const tunnelRouter = router({
  /** Receive the requests (and aborts) for this worker's session */
  serve: publicProcedure
    .input(SessionInput)
    .subscription(async function* ({ ctx, input, signal }): AsyncGenerator<WorkerEvent> {
      await authorizeWorker(ctx, input.session_id)
      const queue: WorkerEvent[] = []
      let wake: (() => void) | undefined
      const notify = () => { wake?.() }
      signal?.addEventListener('abort', notify, { once: true })
      const disconnect = connectWorker(input.session_id, (event) => { queue.push(event); notify() })
      try {
        while (!signal?.aborted) {
          while (queue.length > 0) {
            const event = queue.shift()!
            yield event
            if (event.type === 'close') return
          }
          await new Promise<void>(resolve => { wake = resolve })
          wake = undefined
        }
      } finally {
        signal?.removeEventListener('abort', notify)
        disconnect()
      }
    }),

  /** Respond to a request: a head, then body chunks (base64), then the end */
  respond: publicProcedure
    .input(SessionInput.extend({
      response: z.discriminatedUnion('type', [
        z.object({ type: z.literal('head'), id: z.string(), status: z.number().int(), headers: z.record(z.string(), z.string()) }),
        z.object({ type: z.literal('chunk'), id: z.string(), data: z.string() }),
        z.object({ type: z.literal('end'), id: z.string() }),
        z.object({ type: z.literal('error'), id: z.string(), message: z.string() }),
      ]),
    }))
    .mutation(async ({ ctx, input }) => {
      await authorizeWorker(ctx, input.session_id)
      respond(input.session_id, input.response)
      return null
    }),
})
