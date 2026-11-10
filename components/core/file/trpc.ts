import { z } from 'zod'
import fs from 'fs'
import path from 'path'
import { authedProcedure, router } from '@/app/trpc/init'
import { fetchFile } from './server/fetch'
import { fileFromStream, uploadFile } from './server/upload'
import { abortUpload, completeUpload, startUpload, uploadPartUrls } from './server/direct-upload'
import { waitForVerifiedUpload } from './server/dedupe'
import { TRPCError } from '@trpc/server'
import { NotFoundError } from '@/spec/error'

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

  /** Upload a file from the browser directly to storage, in parts (see ./server/direct-upload.ts) */
  uploadStart: authedProcedure
    .input(z.object({ filename: z.string().min(1).max(1024), size: z.number().int().nonnegative() }))
    .mutation(async ({ ctx, input }) => await startUpload(input, ctx.session.user.id)),

  uploadPartUrls: authedProcedure
    .input(z.object({ token: z.string(), partNumbers: z.array(z.number().int().positive()).min(1).max(100) }))
    .mutation(async ({ ctx, input }) => await uploadPartUrls(input.token, input.partNumbers, ctx.session.user.id)),

  uploadComplete: authedProcedure
    .input(z.object({
      token: z.string(),
      parts: z.array(z.object({ partNumber: z.number().int().positive(), etag: z.string().optional() })),
      sha256: z.string().regex(/^[0-9a-f]{64}$/),
    }))
    .mutation(async ({ ctx, input }) => await completeUpload(input.token, input, ctx.session.user.id)),

  /** Once the upload has been verified, its record (with its `drs://.../<sha256>` url) */
  uploadVerified: authedProcedure
    .input(z.object({ id: z.string() }))
    .subscription(async function* ({ ctx, input, signal }) {
      try {
        yield await waitForVerifiedUpload(input.id, ctx.session.user.id, signal)
      } catch (e) {
        // errors thrown while streaming aren't translated by middleware
        if (e instanceof NotFoundError) throw new TRPCError({ code: 'NOT_FOUND' })
        throw e
      }
    }),

  uploadAbort: authedProcedure
    .input(z.object({ token: z.string() }))
    .mutation(async ({ ctx, input }) => { await abortUpload(input.token, ctx.session.user.id) }),

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
