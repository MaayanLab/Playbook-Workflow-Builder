import { z } from 'zod'
import { TRPCError } from '@trpc/server'
import fpprg from '@/app/fpprg'
import { bus } from '@/app/db'
import { TimeoutError } from '@/spec/error'
import type { Resolved } from '@/core/FPPRG'
import type { StatusUpdate } from '@/spec/metanode'
import { publicProcedure, router } from '../init'

export type ResolvedJSON = ReturnType<Resolved['toJSON']>
export type ResolvedLifecycle = {
  type: 'error',
  error: string,
} | {
  type: 'resolving',
  percent: number,
  status: string,
} | {
  type: 'resolved',
  data: ResolvedJSON
}

export const fpprgRouter = router({
  /** The resolved sequence of steps in a playbook */
  fpl: publicProcedure
    .input(z.string())
    .query(async ({ input: id }) => {
      if (id === 'start') return []
      const fpl = await fpprg.getFPL(id)
      if (!fpl) throw new TRPCError({ code: 'NOT_FOUND', message: 'fpl not found' })
      return fpl.resolve().map(item => item.toJSON())
    }),

  /**
   * Request the resolution of a process, receiving progress updates
   *  until it is resolved.
   */
  resolved: publicProcedure
    .input(z.string())
    .subscription(async function* ({ input: id, signal }): AsyncGenerator<ResolvedLifecycle> {
      const pending: ResolvedLifecycle[] = []
      let outcome: { resolved: Resolved } | { error: unknown } | { done: true } | undefined
      let wake: (() => void) | undefined
      const notify = () => { wake?.() }
      signal?.addEventListener('abort', notify, { once: true })

      // forward partial resolution updates to the client
      const status = { type: 'resolving' as const, percent: 0, status: '' }
      const unsub = bus.listen(`distributed:resolved:${id}:status`, (update: StatusUpdate) => {
        if (update.type === 'progress') status.percent = update.percent
        else if (update.type === 'info') status.status += update.message
        pending.push({ ...status })
        notify()
      })

      // wait for it to be resolved, awaitResolved re-submits the job each time it times out
      ;(async () => {
        while (!signal?.aborted) {
          try {
            return { resolved: await fpprg.awaitResolved(id) }
          } catch (e) {
            if (!TimeoutError.isinstance(e)) return { error: e }
          }
        }
        return { done: true as const }
      })().then(result => { outcome = result; notify() })

      try {
        while (true) {
          while (pending.length > 0) yield pending.shift()!
          if (outcome !== undefined) {
            if ('resolved' in outcome) yield { type: 'resolved', data: outcome.resolved.toJSON() }
            else if ('error' in outcome) {
              console.error(outcome.error)
              yield { type: 'error', error: String(outcome.error) }
            }
            return
          }
          if (signal?.aborted) return
          await new Promise<void>(resolve => { wake = resolve })
          wake = undefined
        }
      } finally {
        signal?.removeEventListener('abort', notify)
        unsub()
      }
    }),
})
