import { z } from 'zod'
import { TRPCError } from '@trpc/server'
import fpprg from '@/app/fpprg'
import { DataC, FPL, IdOrPlaybookMetadataC, IdOrProcessC, Process } from '@/core/FPPRG'
import { decode_complete_process_inputs } from '@/core/engine'
import { getKRG } from '@/app/server-krg'
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

export type ResolvedEvent = { id: string } & ResolvedLifecycle

/**
 * Report the resolution lifecycle of process `id` to `emit` until it's resolved (or fails).
 *  awaitResolved (re-)submits the job each time it times out, so this keeps waiting.
 */
async function watchResolution(id: string, emit: (event: ResolvedLifecycle) => void, signal?: AbortSignal) {
  // forward partial resolution updates
  const status = { type: 'resolving' as const, percent: 0, status: '' }
  const unsub = bus.listen(`distributed:resolved:${id}:status`, (update: StatusUpdate) => {
    if (update.type === 'progress') status.percent = update.percent
    else if (update.type === 'info') status.status += update.message
    emit({ ...status })
  })
  try {
    while (!signal?.aborted) {
      try {
        const resolved = await fpprg.awaitResolved(id)
        emit({ type: 'resolved', data: resolved.toJSON() })
        return
      } catch (e) {
        if (!TimeoutError.isinstance(e)) {
          console.error(e)
          emit({ type: 'error', error: String(e) })
          return
        }
      }
    }
  } finally {
    unsub()
  }
}

export const fpprgRouter = router({
  import: publicProcedure
    .input(z.object({
      data: z.record(z.string(), DataC).optional(),
      workflow: z.array(IdOrProcessC),
      metadata: IdOrPlaybookMetadataC.optional(),
    }))
    .mutation(async ({ input }) => (await fpprg.resolveFPL(input)).id),

  /** Create the first step of a playbook, or extend an existing playbook. */
  extend: publicProcedure
    .input(z.object({ fpl_id: z.string(), process: IdOrProcessC }))
    .mutation(async ({ input }) => {
      const process = await fpprg.resolveProcess(input.process)
      if (input.fpl_id === 'start') return (await fpprg.upsertFPL(new FPL(process))).id
      const oldFpl = await fpprg.getFPL(input.fpl_id)
      if (!oldFpl) throw new TRPCError({ code: 'NOT_FOUND' })
      return (await fpprg.upsertFPL(oldFpl.extend(process))).id
    }),

  rebase: publicProcedure
    .input(z.object({ fpl_id: z.string(), process_id: z.string(), process: IdOrProcessC }))
    .mutation(async ({ input }) => {
      const [oldProcess, oldFpl, newProcess] = await Promise.all([
        fpprg.getProcess(input.process_id),
        fpprg.getFPL(input.fpl_id),
        fpprg.resolveProcess(input.process),
      ])
      if (!oldProcess || !oldFpl) throw new TRPCError({ code: 'NOT_FOUND' })
      const { rebased, head } = oldFpl.rebase(oldProcess, newProcess)
      return { head: (await fpprg.upsertFPL(head)).id, rebased: rebased.id }
    }),

  rebaseMetadata: publicProcedure
    .input(z.object({ fpl_id: z.string(), rebase_id: z.string(), metadata: IdOrPlaybookMetadataC }))
    .mutation(async ({ input }) => {
      const [rebaseFpl, oldFpl, metadata] = await Promise.all([
        fpprg.getFPL(input.rebase_id),
        fpprg.getFPL(input.fpl_id),
        fpprg.resolvePlaybookMetadata(input.metadata),
      ])
      if (!rebaseFpl || !oldFpl) throw new TRPCError({ code: 'NOT_FOUND' })
      const { rebased, head } = oldFpl.rebasePlaybookMetadata(rebaseFpl, metadata)
      return { head: (await fpprg.upsertFPL(head)).id, rebased: rebased.id }
    }),

  recompute: publicProcedure
    .input(z.object({ fpl_id: z.string(), rebase_id: z.string() }))
    .mutation(async ({ input }) => {
      const [rebase, oldFpl] = await Promise.all([fpprg.getFPL(input.rebase_id), fpprg.getFPL(input.fpl_id)])
      if (!rebase || !oldFpl) throw new TRPCError({ code: 'NOT_FOUND' })
      const oldProcess = rebase.process
      if (oldProcess.timestamp !== undefined) {
        const newProcess = await fpprg.upsertProcess(new Process(oldProcess.type, oldProcess.data, oldProcess.inputs))
        const { rebased, head } = oldFpl.rebase(oldProcess, newProcess)
        return { head: (await fpprg.upsertFPL(head)).id, rebased: rebased.id }
      }
      const resolved = await fpprg.getResolved(oldProcess.id)
      if (!resolved) throw new TRPCError({ code: 'NOT_FOUND' })
      await fpprg.deleteResolved(resolved)
      return { head: input.fpl_id, rebased: input.rebase_id }
    }),

  process: publicProcedure
    .input(z.string())
    .query(async ({ input }) => {
      const process = await fpprg.getProcess(input)
      if (!process) throw new TRPCError({ code: 'NOT_FOUND' })
      return process.toJSON()
    }),

  inputs: publicProcedure
    .input(z.string())
    .query(async ({ input }) => {
      const process = await fpprg.getProcess(input)
      if (!process) throw new TRPCError({ code: 'NOT_FOUND' })
      const krg = await getKRG()
      return await decode_complete_process_inputs(krg, process)
    }),

  output: publicProcedure
    .input(z.string())
    .query(async ({ input }) => {
      const process = await fpprg.getProcess(input)
      if (!process) throw new TRPCError({ code: 'NOT_FOUND' })
      return (await process.output())?.toJSON() ?? null
    }),

  deleteOutput: publicProcedure
    .input(z.string())
    .mutation(async ({ input }) => {
      const resolved = await fpprg.getResolved(input)
      if (!resolved) throw new TRPCError({ code: 'NOT_FOUND' })
      await fpprg.deleteResolved(resolved)
      return null
    }),

  fplWithOutput: publicProcedure
    .input(z.string())
    .query(async ({ input }) => {
      const fpl = await fpprg.getFPL(input)
      if (!fpl) throw new TRPCError({ code: 'NOT_FOUND' })
      return await Promise.all(fpl.resolve().map(step => step.toJSONWithOutput()))
    }),


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
  /**
   * Request the resolution of several processes, receiving progress updates
   *  (as `{ id, ...ResolvedLifecycle }`) until all of them are resolved.
   *  A page uses one of these for all of its steps.
   */
  resolvedMany: publicProcedure
    .input(z.array(z.string()).min(1).max(1000))
    .subscription(async function* ({ input, signal }): AsyncGenerator<ResolvedEvent> {
      const ids = [...new Set(input)]
      const queue: ResolvedEvent[] = []
      let wake: (() => void) | undefined
      const notify = () => { wake?.() }
      signal?.addEventListener('abort', notify, { once: true })
      let remaining = ids.length
      for (const id of ids) {
        watchResolution(id, (event) => { queue.push({ id, ...event }); notify() }, signal)
          .finally(() => { remaining -= 1; notify() })
      }
      try {
        while (true) {
          while (queue.length > 0) yield queue.shift()!
          if (remaining === 0 || signal?.aborted) return
          await new Promise<void>(resolve => { wake = resolve })
          wake = undefined
        }
      } finally {
        signal?.removeEventListener('abort', notify)
      }
    }),
})
