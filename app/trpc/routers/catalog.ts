import { TRPCError } from '@trpc/server'
import { z } from 'zod'
import krg from '@/app/krg'
import { metanodeType } from '@/app/fragments/components/listing'
import * as dict from '@/utils/dict'
import { publicProcedure, router } from '../init'

function summary(metanode: ReturnType<typeof krg.getDataNodes>[number] | ReturnType<typeof krg.getProcessNodes>[number]) {
  const { spec: id, meta: { icon: _, ...meta }, kind } = metanode
  return { id, kind, type: metanodeType(metanode), meta }
}

export const catalogRouter = router({
  resolve: publicProcedure
    .input(z.object({
      spec: z.string(),
      data: z.unknown().optional(),
      inputs: z.record(z.string(), z.unknown()),
    }))
    .mutation(async ({ input }) => {
      const process = krg.getProcessNode(input.spec)
      const output = await process.resolve({ data: input.data, inputs: input.inputs, notify: () => {} })
      return process.output.codec.encode(output)
    }),

  list: publicProcedure
    .query(() => [...krg.getDataNodes(), ...krg.getProcessNodes()].map(summary)),

  get: publicProcedure
    .input(z.string())
    .query(({ input: spec }) => {
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
    .query(({ input }) => krg.getNextProcess(input).map(summary)),

  previous: publicProcedure
    .input(z.string())
    .query(({ input }) => krg.getPrevProcess(input).map(summary)),
})
