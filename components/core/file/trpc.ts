import { z } from 'zod'
import fs from 'fs'
import path from 'path'
import { authedProcedure, router } from '@/app/trpc/init'
import { fetchFile } from './server/fetch'
import { fileFromStream, uploadFile } from './server/upload'
import { abortUpload, completeUpload, startUpload, uploadPartUrls } from './server/direct-upload'
import { uploadStatus, waitForVerifiedUpload } from './server/dedupe'
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

  /**
   * Upload a file directly to storage, in parts (see ./server/direct-upload.ts), used by the
   *  browser & documented for API users (docs/user/api.md, as REST endpoints under /api/v1/uploads)
   */
  uploadStart: authedProcedure
    .meta({ openapi: {
      method: 'POST', path: '/v1/uploads', tags: ['Uploads'], protect: true,
      summary: 'Start uploading a file',
      description: 'Returns a token for the upload & how to split the file: `partCount` parts of `partSize` bytes (the last may be smaller).',
    } })
    .input(z.object({ filename: z.string().min(1).max(1024), size: z.number().int().nonnegative() }))
    .output(z.object({ token: z.string(), partSize: z.number(), partCount: z.number(), etagsRequired: z.boolean() }))
    .mutation(async ({ ctx, input }) => await startUpload(input, ctx.session.user.id)),

  uploadPartUrls: authedProcedure
    .meta({ openapi: {
      method: 'POST', path: '/v1/uploads/parts', tags: ['Uploads'], protect: true,
      summary: 'Get urls to upload parts to',
      description: 'PUT each part (its exact bytes) to its `url`, with its `headers`, and keep the `ETag` response header. At most 100 parts per request, urls expire after an hour.',
    } })
    .input(z.object({ token: z.string(), partNumbers: z.array(z.number().int().positive()).min(1).max(100) }))
    .output(z.array(z.object({ partNumber: z.number(), url: z.string(), headers: z.record(z.string(), z.string()), proxy: z.boolean() })))
    .mutation(async ({ ctx, input }) => await uploadPartUrls(input.token, input.partNumbers, ctx.session.user.id)),

  uploadComplete: authedProcedure
    .meta({ openapi: {
      method: 'POST', path: '/v1/uploads/complete', tags: ['Uploads'], protect: true,
      summary: 'Complete an upload',
      description: "Assembles the file from its parts. It's then verified & deduplicated in the background, see `GET /v1/uploads/{id}`.",
    } })
    .input(z.object({
      token: z.string(),
      parts: z.array(z.object({ partNumber: z.number().int().positive(), etag: z.string().optional() })),
      sha256: z.string().regex(/^[0-9a-f]{64}$/),
    }))
    .output(z.object({ id: z.string() }))
    .mutation(async ({ ctx, input }) => await completeUpload(input.token, input, ctx.session.user.id)),

  uploadStatus: authedProcedure
    .meta({ openapi: {
      method: 'GET', path: '/v1/uploads/{id}', tags: ['Uploads'], protect: true,
      summary: "An upload's status",
      description: "`verifying` until the file is verified, then `verified` with the file's record: its `url` (`drs://.../<sha256>`) refers to the file in workflows.",
    } })
    .input(z.object({ id: z.string() }))
    .output(z.object({
      status: z.enum(['verifying', 'verified']),
      id: z.string(), url: z.string().optional(), filename: z.string().optional(), sha256: z.string().optional(), size: z.number().optional(),
    }))
    .query(async ({ ctx, input }) => await uploadStatus(input.id, ctx.session.user.id)),

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
    .meta({ openapi: { method: 'POST', path: '/v1/uploads/abort', tags: ['Uploads'], protect: true, summary: 'Abort an upload' } })
    .input(z.object({ token: z.string() }))
    .output(z.object({}))
    .mutation(async ({ ctx, input }) => { await abortUpload(input.token, ctx.session.user.id); return {} }),

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
