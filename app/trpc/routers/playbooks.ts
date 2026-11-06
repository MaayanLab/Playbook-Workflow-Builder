import { z } from 'zod'
import { sql } from 'kysely'
import { TRPCError } from '@trpc/server'
import db from '@/app/db'
import krg from '@/app/krg'
import fpprg from '@/app/fpprg'
import { uuid, SYSTEM_USER_ID, type UserPlaybook } from '@/db'
import * as dict from '@/utils/dict'
import { IdOrPlaybookMetadataC, IdOrCellMetadataC } from '@/core/FPPRG'
import { tsvector, tsvector_intersect } from '@/utils/tsvector'
import cache from '@/utils/global_cache'
import publicPlaybooks from '@/app/public/playbooksDemo'
import { authedProcedure, publicProcedure, router } from '../init'

const playbook_tsvectors = cache('playbook_tsvectors', () => {
  const playbook_tsvectors: Record<string, Set<string>> = {}
  publicPlaybooks.forEach(playbook => {
    playbook.workflow.workflow[playbook.workflow.workflow.length-1].id = playbook.id
    playbook_tsvectors[playbook.id] = tsvector([
      playbook.label,
      playbook.description,
      ...playbook.inputs.flatMap(input => [
        input.meta.label,
        input.meta.description,
      ]),
      ...playbook.outputs.flatMap(output => [
        output.meta.label,
        output.meta.description,
      ]),
      ...playbook.dataSources,
    ].join(' '))
  })
  return playbook_tsvectors
})

const notFound = () => new TRPCError({ code: 'NOT_FOUND' })

export const playbooksRouter = router({
  /** Curated playbooks shipped with the app */
  public: publicProcedure
    .input(z.object({
      search: z.string().optional(),
      inputs: z.string().optional(),
      outputs: z.string().optional(),
      dataSources: z.string().optional(),
      skip: z.number().optional(),
      limit: z.number().optional(),
    }))
    .query(async ({ input }) => {
      const search = input.search
      const inputs = input.inputs ? input.inputs.split(', ') : undefined
      const outputs = input.outputs ? input.outputs.split(', ') : undefined
      const dataSources = input.dataSources ? input.dataSources.split(', ') : undefined
      const skip = input.skip ?? 0
      const limit = input.limit ?? 50
      let playbooks = publicPlaybooks
      if (inputs) playbooks = playbooks.filter(playbook => !inputs.some(spec => !playbook.inputs.map(t=>t.spec as string).includes(spec)))
      if (outputs) playbooks = playbooks.filter(playbook => !outputs.some(spec => !playbook.outputs.map(t=>t.spec as string).includes(spec)))
      if (dataSources) playbooks = playbooks.filter(playbook => !dataSources.some(dataSource => !playbook.dataSources.includes(dataSource)))
      if (search) {
        const search_tsvector = tsvector(search)
        const search_scores: Record<string, number> = {}
        playbooks.forEach(playbook => {
          search_scores[playbook.id] = tsvector_intersect(playbook_tsvectors[playbook.id], search_tsvector).size
        })
        playbooks = playbooks.filter(playbook => search_scores[playbook.id] > 0)
        playbooks.sort((a, b) => search_scores[b.id] - search_scores[a.id])
      }
      playbooks = playbooks.slice(skip, limit)
      const userPlaybooks = playbooks.length === 0 ? [] : await db.selectFrom('user_playbook')
        .selectAll()
        .where('playbook', 'in', playbooks.map(({ id }) => id))
        .orderBy('created', 'desc')
        .execute()
      const userPlaybookLookup: Record<string, UserPlaybook | undefined> = dict.init(userPlaybooks.map((playbook) => ({ key: playbook.playbook, value: playbook })))
      const results = await Promise.all(playbooks.map(async ({ inputs, outputs, dataSources, ...playbook }) => {
        if (userPlaybookLookup[playbook.id] === undefined) {
          try {
            const fpl = await fpprg.resolveFPL(playbook.workflow as any)
            fpl.id = playbook.id
            await fpprg.upsertFPL(fpl)
            // registered as public so that views are counted (see `get`)
            userPlaybookLookup[playbook.id] = await db.insertInto('user_playbook')
              .values({
                id: uuid(),
                user: SYSTEM_USER_ID,
                playbook: playbook.id,
                title: playbook.label,
                description: playbook.description,
                public: true,
              })
              .returningAll()
              .executeTakeFirstOrThrow()
          } catch (e) {
            console.error(e)
          }
        }
        return {
          ...playbook,
          inputs: inputs.map(t => t.spec as string).join(', '),
          outputs: outputs.map(t => t.spec as string).join(', '),
          dataSources: dataSources.join(', '),
          clicks: userPlaybookLookup[playbook.id]?.clicks || 0,
          disabled: userPlaybookLookup[playbook.id] === undefined,
        }
      }))
      if (!search) results.sort((a, b) => b.clicks - a.clicks)
      return results
    }),

  /** Playbooks published by users */
  community: publicProcedure
    .input(z.object({
      search: z.string().optional(),
      inputs: z.string().optional(),
      outputs: z.string().optional(),
      skip: z.number().optional(),
      limit: z.number().optional(),
    }))
    .query(async ({ input }) => {
      const search = input.search
      const inputs = input.inputs ? input.inputs.split(', ') : undefined
      const outputs = input.outputs ? input.outputs.split(', ') : undefined
      // TODO: filter in DB
      let playbooks = await db.selectFrom('user_playbook')
        .selectAll()
        .where('public', '=', true)
        .where('user', '!=', SYSTEM_USER_ID)
        .orderBy('clicks', 'desc')
        .offset(input.skip ?? 0)
        .limit(input.limit ?? 50)
        .execute()
      if (inputs) playbooks = playbooks.filter(playbook => !inputs.some(spec => !(playbook.inputs||'').split(', ').includes(spec)))
      if (outputs) playbooks = playbooks.filter(playbook => !outputs.some(spec => !(playbook.outputs||'').split(', ').includes(spec)))
      if (search) {
        const search_tsvector = tsvector(search)
        const search_scores = dict.init(playbooks.map(playbook => {
          const playbook_tsvector = tsvector([
            playbook.title || '',
            playbook.description || '',
          ].join(' '))
          return { key: playbook.id, value: tsvector_intersect(search_tsvector, playbook_tsvector).size }
        }))
        playbooks.sort((a, b) => search_scores[b.id] - search_scores[a.id])
      }
      return playbooks
    }),

  /** The current user's saved playbooks */
  mine: authedProcedure
    .query(async ({ ctx }) => {
      return await db.selectFrom('user_playbook')
        .selectAll()
        .where('user', '=', ctx.session.user.id)
        .execute()
    }),

  /** Record a view of the playbook & get the current user's sharing status for it */
  get: publicProcedure
    .input(z.object({ id: z.string() }))
    .query(async ({ ctx, input }) => {
      const fpl = await fpprg.getFPL(input.id)
      if (fpl === undefined) throw notFound()
      // TODO: move this logic into a db trigger
      // TODO: address playbook republishing
      const publicUserPlaybook = await db.selectFrom('user_playbook')
        .select('id')
        .where('playbook', '=', input.id)
        .where('public', '=', true)
        .orderBy('created', 'asc')
        .limit(1)
        .executeTakeFirst()
      if (publicUserPlaybook) {
        await db.updateTable('user_playbook')
          .set({ clicks: sql`clicks + 1` })
          .where('id', '=', publicUserPlaybook.id)
          .execute()
      }
      const session = await ctx.getSession()
      if (!session?.user) return null
      const userPlaybook = await db.selectFrom('user_playbook')
        .select('public')
        .where('user', '=', session.user.id)
        .where('playbook', '=', input.id)
        .executeTakeFirst()
      if (!userPlaybook) return null
      return { public: userPlaybook.public }
    }),

  /** Save the playbook (with updated metadata) to the current user's account */
  update: authedProcedure
    .input(z.object({
      id: z.string(),
      user_playbook: z.object({ public: z.boolean() }),
      playbook_metadata: IdOrPlaybookMetadataC,
      cell_metadata: z.record(z.string(), IdOrCellMetadataC),
    }))
    .mutation(async ({ ctx, input }) => {
      let fpl = await fpprg.getFPL(input.id)
      if (fpl === undefined) throw notFound()
      fpl = fpl.rebaseCellMetadata(
        dict.init(
          await Promise.all(dict.items(input.cell_metadata)
            .map(async ({ key, value }) => ({ key, value: await fpprg.resolveCellMetadata(value) })))
        )
      )
      if (fpl === undefined) throw notFound()
      fpl = fpl.rebasePlaybookMetadata(fpl, await fpprg.resolvePlaybookMetadata(input.playbook_metadata)).head
      fpl = await fpprg.upsertFPL(fpl)
      if (fpl === undefined) throw notFound()
      fpl = await fpl.rebaseInTime()
      if (fpl === undefined) throw notFound()
      fpl = await fpprg.upsertFPL(fpl)
      const playbookInputs = fpl.resolve().map(cell => krg.getProcessNode(cell.process.type)).filter(node => 'prompt' in node && dict.isEmpty(node.inputs)).map(node => node.output.spec).join(', ')
      const playbookOutputs = fpl.resolve().filter(cell => cell.cell_metadata?.data_visible === true).map(cell => krg.getProcessNode(cell.process.type).output.spec).join(', ')
      const record = {
        title: fpl.playbook_metadata?.title ?? null,
        description: fpl.playbook_metadata?.description ?? null,
        public: input.user_playbook.public,
        inputs: playbookInputs,
        outputs: playbookOutputs,
      }
      const user = ctx.session.user.id
      const playbook = fpl.id
      await db.transaction().execute(async (trx) => {
        // the user's entry for the previous version is replaced
        await trx.deleteFrom('user_playbook')
          .where('user', '=', user)
          .where('playbook', '=', input.id)
          .execute()
        const existing = await trx.updateTable('user_playbook')
          .set(record)
          .where('user', '=', user)
          .where('playbook', '=', playbook)
          .executeTakeFirst()
        if (Number(existing.numUpdatedRows) === 0) {
          await trx.insertInto('user_playbook')
            .values({ id: uuid(), user, playbook, ...record })
            .execute()
        }
      })
      return fpl.id
    }),

  publish: authedProcedure
    .input(z.object({ id: z.string(), public: z.boolean() }))
    .mutation(async ({ ctx, input }) => {
      const userPlaybook = await db.updateTable('user_playbook')
        .set({ public: input.public })
        .where('user', '=', ctx.session.user.id)
        .where('playbook', '=', input.id)
        .returningAll()
        .executeTakeFirst()
      return userPlaybook ?? null
    }),

  delete: authedProcedure
    .input(z.object({ id: z.string() }))
    .mutation(async ({ ctx, input }) => {
      const fpl = await fpprg.getFPL(input.id)
      if (fpl === undefined) throw notFound()
      const userPlaybook = await db.deleteFrom('user_playbook')
        .where('user', '=', ctx.session.user.id)
        .where('playbook', '=', input.id)
        .returningAll()
        .executeTakeFirst()
      return userPlaybook ?? null
    }),
})
