/**
 * Uploads from the browser are stored under a random key (see ./direct-upload.ts), as the
 *  checksum the browser reports can't be trusted to decide where they're stored. Afterwards, the
 *  workers verify & deduplicate them:
 *
 * 1. The stored object is read to compute its actual sha256 (correcting the upload's, if needed).
 * 2. The upload is pointed at `storage://<sha256>`, copying the object there (in storage) unless
 *    an identical file is already stored.
 * 3. The random key copy is removed later, once urls to it that were already given out
 *    (e.g. presigned urls) have expired, provided nothing refers to it anymore.
 */
import crypto from 'crypto'
import db, { bus } from '@/app/db'
import { NotFoundError } from '@/spec/error'
import storage, { storageKey } from './storage'
import { drsUrl } from './drs'

const DEDUPE_QUEUE = 'upload-dedupe-queue'
const DELETE_QUEUE = 'storage-delete-queue'
/** announces an upload's verification: `{ error?: string }` */
const uploadEvent = (uploadId: string) => `distributed:upload:${uploadId}`
/** longer than presigned urls are valid (an hour) */
const deleteDelaySeconds = () => +(process.env.STORAGE_DELETE_DELAY_SECONDS || 2 * 60 * 60)

/** Verify & deduplicate upload `uploadId` (an `upload` row) in the background */
export async function queueDedupe(uploadId: string) {
  // the workers start with the fpprg (where configured), without a database the work is
  //  dispatched in process & would otherwise be dropped if they hadn't been started yet
  await import('@/app/fpprg')
  await bus.send(DEDUPE_QUEUE, { id: uploadId })
}

const CHUNK_SIZE = 64 * 1024 * 1024
const CHUNK_ATTEMPTS = 3

/** Read the (inclusive) byte range of an object, retrying if reading it fails */
async function readChunk(key: string, start: number, end: number) {
  for (let attempt = 1; ; attempt++) {
    try {
      const chunks: Buffer[] = []
      for await (const chunk of await storage().read(key, { start, end })) chunks.push(Buffer.from(chunk))
      const buffer = Buffer.concat(chunks)
      if (buffer.length !== end - start + 1) throw new Error(`read ${buffer.length} bytes of ${start}-${end}`)
      return buffer
    } catch (e) {
      if (attempt >= CHUNK_ATTEMPTS) throw e
      await new Promise(resolve => setTimeout(resolve, 1000 * attempt))
    }
  }
}

/**
 * The sha256 & size of a stored object, read in chunks (so that a failure only repeats a chunk),
 *  the next chunk is read while the current one is hashed.
 */
export async function contentStats(key: string) {
  const size = (await storage().stat(key))?.size
  if (size === undefined) throw new NotFoundError()
  const hash = crypto.createHash('sha256')
  const read = (start: number) => {
    const chunk = readChunk(key, start, Math.min(start + CHUNK_SIZE, size) - 1)
    // (it's awaited later, a failure shouldn't be reported as unhandled before then)
    chunk.catch(() => {})
    return chunk
  }
  let next = size > 0 ? read(0) : undefined
  for (let start = 0; start < size; start += CHUNK_SIZE) {
    const chunk = await next!
    next = start + CHUNK_SIZE < size ? read(start + CHUNK_SIZE) : undefined
    hash.update(chunk)
  }
  return { sha256: hash.digest('hex'), size }
}

async function dedupe(uploadId: string) {
  const upload = await db.selectFrom('upload').selectAll().where('id', '=', uploadId).executeTakeFirst()
  if (!upload || !upload.url.startsWith('storage://')) return
  const key = storageKey(upload.url)
  const { sha256, size } = await contentStats(key)
  // already content addressed
  if (key === sha256) return
  if (sha256 !== upload.sha256 || size !== upload.size) {
    console.warn(`upload ${uploadId} was recorded as ${upload.sha256} (${upload.size} bytes) but is ${sha256} (${size} bytes), corrected`)
  }
  const backend = storage()
  if (!(await backend.stat(sha256))) await backend.copy(key, sha256)
  await db.updateTable('upload').set({ url: `storage://${sha256}`, sha256, size }).where('id', '=', uploadId).execute()
  await bus.send(DELETE_QUEUE, { id: key, delaySeconds: deleteDelaySeconds() })
  console.debug(`upload ${uploadId} deduplicated as storage://${sha256}`)
}

/** Remove an object from storage, unless an upload refers to it (again) */
async function deleteUnreferenced(key: string) {
  const referenced = await db.selectFrom('upload').select('id').where('url', '=', `storage://${key}`).executeTakeFirst()
  if (referenced) return
  await storage().delete(key)
  console.debug(`removed storage://${key}`)
}

export function start_dedupe_workers(concurrency: number) {
  bus.work(DEDUPE_QUEUE, { concurrency }, async (job) => {
    try {
      await dedupe(job.data.id)
      bus.notify(uploadEvent(job.data.id), {})
    } catch (e) {
      console.error(`deduplicating upload ${job.data.id} failed`, e)
      bus.notify(uploadEvent(job.data.id), { error: e instanceof Error ? e.message : String(e) })
    }
  })
  bus.work(DELETE_QUEUE, { concurrency }, async (job) => {
    try { await deleteUnreferenced(job.data.id) }
    catch (e) { console.error(`removing storage://${job.data.id} failed`, e) }
  })
}

export type VerifiedUpload = { id: string, url: string, filename: string, sha256: string, size: number }

/**
 * Wait for user upload `userUploadId` (of `user`) to be verified, giving its record
 *  (with its content addressed url) once it has been.
 */
export async function waitForVerifiedUpload(userUploadId: string, user: string, signal?: AbortSignal): Promise<VerifiedUpload> {
  const row = await db.selectFrom('user_upload').select(['upload', 'filename']).where('id', '=', userUploadId).where('user', '=', user).executeTakeFirst()
  if (!row) throw new NotFoundError()
  // listen before checking, so that no notification can be missed in between
  let notified: { error?: string } | undefined
  let wake: (() => void) | undefined
  const unlisten = bus.listen(uploadEvent(row.upload), (data: { error?: string }) => { notified = data; wake?.() })
  const onAbort = () => wake?.()
  signal?.addEventListener('abort', onAbort, { once: true })
  try {
    while (!signal?.aborted) {
      const upload = await db.selectFrom('upload').select(['url', 'sha256', 'size']).where('id', '=', row.upload).executeTakeFirstOrThrow()
      if (upload.url === `storage://${upload.sha256}`) {
        return { id: userUploadId, url: drsUrl(upload.sha256), filename: row.filename, sha256: upload.sha256, size: upload.size }
      }
      if (notified?.error) throw new Error(`The upload couldn't be verified: ${notified.error}`)
      // re-check periodically in case a notification was lost
      await new Promise<void>(resolve => {
        const timeout = setTimeout(resolve, 15_000)
        wake = () => { clearTimeout(timeout); resolve() }
      })
      wake = undefined
    }
    throw new Error('aborted')
  } finally {
    unlisten()
    signal?.removeEventListener('abort', onAbort)
  }
}

/** Whether user upload `userUploadId` (of `user`) has been verified (& its record if so) */
export async function uploadStatus(userUploadId: string, user: string) {
  const row = await db.selectFrom('user_upload')
    .innerJoin('upload', 'upload.id', 'user_upload.upload')
    .select(['user_upload.filename', 'upload.url', 'upload.sha256', 'upload.size'])
    .where('user_upload.id', '=', userUploadId)
    .where('user_upload.user', '=', user)
    .executeTakeFirst()
  if (!row) throw new NotFoundError()
  if (row.url !== `storage://${row.sha256}`) return { status: 'verifying' as const, id: userUploadId }
  return { status: 'verified' as const, id: userUploadId, url: drsUrl(row.sha256), filename: row.filename, sha256: row.sha256, size: row.size }
}
