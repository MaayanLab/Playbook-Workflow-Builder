/**
 * Exports produce files which are kept (as the user's uploads), so they require signing in,
 *  and users can only see their own exports.
 */
import { z } from 'zod'
import { TRPCError } from '@trpc/server'
import fpprg from '@/app/fpprg'
import { ExportFormatC, ExportMetadataC, exportAuthor, getExport, requestExport, watchExport, type ExportState } from '@/app/exports'
import { NotFoundError } from '@/spec/error'
import { authedProcedure, router } from '../init'

export const exportsRouter = router({
  /** Queue an export of a playbook, returning its id */
  create: authedProcedure
    .input(z.object({
      fpl_id: z.string(),
      format: ExportFormatC,
      metadata: ExportMetadataC.optional(),
    }))
    .mutation(async ({ ctx, input }) => {
      const fpl = await fpprg.getFPL(input.fpl_id)
      if (!fpl) throw new NotFoundError()
      const user = ctx.session.user.id
      // (only the LaTeX & PDF exports name an author)
      const author = input.format.startsWith('tex-') ? await exportAuthor(user) : undefined
      const id = await requestExport({
        fpl_id: fpl.id,
        format: input.format,
        params: { metadata: input.metadata, author },
        user,
      })
      return { id }
    }),

  get: authedProcedure
    .input(z.object({ id: z.string() }))
    .query(async ({ ctx, input }): Promise<ExportState> => {
      // other users' exports are indistinguishable from missing ones
      const state = await getExport(input.id, ctx.session.user.id)
      if (!state) throw new NotFoundError()
      return state
    }),

  /** The state of an export each time it changes, ending once it's complete or has failed */
  watch: authedProcedure
    .input(z.object({ id: z.string() }))
    .subscription(async function* ({ ctx, input, signal }) {
      // errors thrown while streaming aren't translated by middleware, so this must be a TRPCError
      if (!await getExport(input.id, ctx.session.user.id)) throw new TRPCError({ code: 'NOT_FOUND' })
      yield* watchExport(input.id, ctx.session.user.id, signal)
    }),
})
