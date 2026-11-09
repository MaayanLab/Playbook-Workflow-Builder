import { TRPCError } from '@trpc/server'
import { z } from 'zod'
import { getComponentKRG } from '@/app/server-krg'
import db from '@/app/db'
import { metanodeType } from '@/app/metanode-type'
import * as dict from '@/utils/dict'
import { publicProcedure, router } from '../init'

type CatalogKRG = ReturnType<typeof getComponentKRG>
type CatalogMetaNode = ReturnType<CatalogKRG['getDataNodes']>[number]
  | ReturnType<CatalogKRG['getProcessNodes']>[number]

function summary(metanode: CatalogMetaNode) {
  const { spec: id, meta: { icon: _, ...meta }, kind } = metanode
  return { id, kind, type: metanodeType(getComponentKRG(), metanode), meta }
}

/**
 * The components catalog (user suggestions are not components, see `suggestions`)
 */
export const catalogRouter = router({
  /** Every user suggestion, the browser adds them to its KRG (see @/app/extensions/suggestions) */
  suggestions: publicProcedure
    .query(async () => await db.selectFrom('suggestion')
      .select(['id', 'name', 'inputs', 'output', 'description', 'user', 'created'])
      .orderBy('created')
      .execute()),

  resolve: publicProcedure
    .input(z.object({
      spec: z.string(),
      data: z.unknown().optional(),
      inputs: z.record(z.string(), z.unknown()),
    }))
    .mutation(async ({ input }) => {
      const krg = getComponentKRG()
      const process = krg.getProcessNode(input.spec)
      if (!process) throw new TRPCError({ code: 'NOT_FOUND' })
      const output = await process.resolve({ data: input.data, inputs: input.inputs, notify: () => {} })
      return process.output.codec.encode(output)
    }),

  list: publicProcedure
    .query(async () => {
      const krg = getComponentKRG()
      return [...krg.getDataNodes(), ...krg.getProcessNodes()].map(summary)
    }),

  get: publicProcedure
    .input(z.string())
    .query(async ({ input: spec }) => {
      const krg = getComponentKRG()
      const data = krg.getDataNode(spec)
      if (data) return summary(data)
      const process = krg.getProcessNode(spec)
      if (!process) throw new TRPCError({ code: 'NOT_FOUND' })
      return {
        ...summary(process),
        inputs: dict.init(dict.items(process.inputs).map(({ key, value }) => ({
          key,
          value: Array.isArray(value) ? value.map(input => ({ id: input.spec })) : { id: value.spec },
        }))),
        output: { id: process.output.spec },
        story: process.story({}),
      }
    }),

  next: publicProcedure
    .input(z.string())
    .query(({ input }) => getComponentKRG().getNextProcess(input).map(summary)),

  previous: publicProcedure
    .input(z.string())
    .query(({ input }) => getComponentKRG().getPrevProcess(input).map(summary)),
})
