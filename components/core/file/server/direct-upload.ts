/**
 * Uploads from the browser directly to storage, in parts:
 *
 * 1. `start` begins a multipart upload & returns a signed token describing it (key, size, owner),
 *    so no server side state is needed in between calls.
 * 2. `partUrls` gives out urls for parts: presigned urls to storage itself where the backend has
 *    them (S3, CAVATICA), or this app's `/api/file/upload/part` (the filesystem, or as a fallback
 *    when a direct PUT fails, e.g. because of CORS).
 * 3. `complete` assembles the object & records the upload.
 * 4. `uploadVerified` waits for the upload to be verified, giving its (content addressed) url.
 *
 * Files uploaded this way are stored under a random key rather than their sha256: the browser
 *  computes the checksum as it uploads, which can't be verified before the data is stored, so
 *  it mustn't determine where (i.e. as which content) it's stored. The size is verified, the
 *  workers then verify the checksum & deduplicate the file (see ./dedupe.ts), after which the
 *  browser is given the file's url: `drs://.../<sha256>` (identical files are identical inputs).
 */
import crypto from 'crypto'
import { z } from 'zod'
import type { Readable } from 'stream'
import db from '@/app/db'
import { insertUserUpload } from '@/db/queries'
import { ResponseCodedError, UnauthorizedError } from '@/spec/error'
import storage, { type MultipartPart } from './storage'
import { queueDedupe } from './dedupe'

const TokenC = z.object({
  key: z.string(),
  uploadId: z.string(),
  size: z.number(),
  partSize: z.number(),
  filename: z.string(),
  user: z.string(),
  exp: z.number(),
})
type Token = z.infer<typeof TokenC>

/** uploads must complete within a day */
const TOKEN_TTL_MS = 24 * 60 * 60 * 1000

let processSecret: string | undefined
function secret() {
  const configured = process.env.BETTER_AUTH_SECRET || process.env.NEXTAUTH_SECRET || process.env.WORKER_AUTH_TOKEN
  if (configured) return configured
  // without a configured secret tokens are only valid within this process
  return processSecret ??= crypto.randomBytes(32).toString('hex')
}

function sign(token: Token) {
  const payload = Buffer.from(JSON.stringify(token)).toString('base64url')
  const mac = crypto.createHmac('sha256', secret()).update(payload).digest('base64url')
  return `${payload}.${mac}`
}

export function verifyUploadToken(signed: string, user?: string) {
  const [payload, mac] = signed.split('.')
  if (!payload || !mac) throw new UnauthorizedError()
  const expected = crypto.createHmac('sha256', secret()).update(payload).digest('base64url')
  if (mac.length !== expected.length || !crypto.timingSafeEqual(Buffer.from(mac), Buffer.from(expected))) throw new UnauthorizedError()
  const token = TokenC.parse(JSON.parse(Buffer.from(payload, 'base64url').toString()))
  if (token.exp < Date.now()) throw new ResponseCodedError(410, 'This upload has expired')
  if (user !== undefined && token.user !== user) throw new UnauthorizedError()
  return token
}

const partCount = (token: Pick<Token, 'size' | 'partSize'>) => Math.max(1, Math.ceil(token.size / token.partSize))
/** the exact size of a part (all are partSize except the last) */
function partLength(token: Token, partNumber: number) {
  const count = partCount(token)
  if (!Number.isInteger(partNumber) || partNumber < 1 || partNumber > count) throw new ResponseCodedError(400, `Invalid part number ${partNumber}`)
  return partNumber < count ? token.partSize : token.size - token.partSize * (count - 1)
}

export async function startUpload({ filename, size }: { filename: string, size: number }, user: string) {
  const key = crypto.randomUUID()
  const { uploadId, partSize } = await storage().createMultipart(key, size)
  const token: Token = { key, uploadId, size, partSize, filename, user, exp: Date.now() + TOKEN_TTL_MS }
  return {
    token: sign(token), partSize, partCount: partCount(token),
    // CAVATICA must be told every part's etag, which the browser can only read if storage allows it
    etagsRequired: storage().type === 'cavatica',
  }
}

export async function uploadPartUrls(signed: string, partNumbers: number[], user: string) {
  const token = verifyUploadToken(signed, user)
  return await Promise.all(partNumbers.map(async (partNumber) => {
    partLength(token, partNumber)
    const direct = await storage().partUrl(token.key, token.uploadId, partNumber)
    return direct ? { partNumber, ...direct, proxy: false as const } : { partNumber, proxy: true as const }
  }))
}

/** A part uploaded through this app (`/api/file/upload/part`) */
export async function uploadPartThroughApp(signed: string, partNumber: number, body: Readable, contentLength: number | undefined) {
  const token = verifyUploadToken(signed)
  const size = partLength(token, partNumber)
  if (contentLength !== undefined && contentLength !== size) throw new ResponseCodedError(400, `Part ${partNumber} must be ${size} bytes`)
  return await storage().uploadPart(token.key, token.uploadId, partNumber, body, size)
}

export async function completeUpload(signed: string, { parts, sha256 }: { parts: MultipartPart[], sha256: string }, user: string) {
  const token = verifyUploadToken(signed, user)
  const count = partCount(token)
  const numbers = new Set(parts.map(part => part.partNumber))
  if (numbers.size !== count || parts.some(part => part.partNumber < 1 || part.partNumber > count)) {
    throw new ResponseCodedError(400, `Expected parts 1-${count}`)
  }
  const backend = storage()
  await backend.completeMultipart(token.key, token.uploadId, parts)
  const stat = await backend.stat(token.key)
  if (stat?.size !== token.size) {
    await backend.delete(token.key).catch(() => {})
    throw new ResponseCodedError(400, `The uploaded file is ${stat?.size ?? 0} bytes, expected ${token.size}`)
  }
  const url = `storage://${token.key}`
  const user_upload = await insertUserUpload(db, { url, sha256, size: token.size, filename: token.filename }, token.user)
  // the workers verify the checksum & deduplicate the file (see ./dedupe.ts)
  const { upload } = await db.selectFrom('user_upload').select('upload').where('id', '=', user_upload.id).executeTakeFirstOrThrow()
  await queueDedupe(upload)
  // the file's (content addressed) url is given out once it's verified (`uploadVerified`)
  return { id: user_upload.id }
}

export async function abortUpload(signed: string, user: string) {
  const token = verifyUploadToken(signed, user)
  await storage().abortMultipart(token.key, token.uploadId)
}
