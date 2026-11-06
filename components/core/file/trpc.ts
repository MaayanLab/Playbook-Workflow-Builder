import { z } from 'zod'
import fs from 'fs'
import path from 'path'
import { authedProcedure, router } from '@/app/trpc/init'
import { fetchFile } from './server/fetch'
import { fileFromStream, uploadFile } from './server/upload'

const example = z.object({
  asset: z.string()
    .regex(/^[A-Za-z0-9._/-]+\/assets\/[A-Za-z0-9._/-]+$/)
    .refine(asset => !asset.split('/').some(part => part === '.' || part === '..'), 'Invalid component asset path'),
  filename: z.string(),
})

/** Typed API owned by the core file component. */
const fileRouter = router({
  fetch: authedProcedure
    .input(z.object({ url: z.string() }))
    .mutation(async ({ ctx, input }) => await fetchFile(input, ctx.session)),

  register: authedProcedure
    .input(example)
    .mutation(async ({ ctx, input }) => {
      const componentsRoot = path.resolve(process.env.APP_ROOT as string, 'components')
      const asset = path.resolve(componentsRoot, input.asset)
      if (path.relative(componentsRoot, asset).startsWith('..')) throw new Error('Invalid component asset path')
      const temporaryFile = await fileFromStream(fs.createReadStream(asset), input.filename)
      return await uploadFile(temporaryFile, ctx.session)
    }),
})

export default fileRouter
