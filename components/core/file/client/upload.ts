/**
 * Upload a file from the browser directly to storage, in parts (see ../server/direct-upload.ts).
 *  The file's sha256 is computed as it's read for uploading, the server then verifies it before
 *  giving the file's (content addressed) url.
 */
import { createSHA256 } from 'hash-wasm'
import type { TRPCClient } from '@trpc/client'
import type { AppRouter } from '@/app/trpc/router'
import type { UploadFileResponse } from '../server/upload'

/** parts uploaded at once */
const CONCURRENCY = 4
/** attempts for a part sent through the app */
const ATTEMPTS = 3

export type UploadClient = TRPCClient<AppRouter>

/** A network level failure (status 0), e.g. storage not allowing this origin (CORS) */
class NetworkError extends Error {}

function put(url: string, headers: Record<string, string>, body: ArrayBuffer, onProgress: (loaded: number) => void) {
  return new Promise<string | null>((resolve, reject) => {
    const xhr = new XMLHttpRequest()
    xhr.open('PUT', url)
    for (const [key, value] of Object.entries(headers)) xhr.setRequestHeader(key, value)
    xhr.upload.onprogress = (evt) => onProgress(evt.loaded)
    xhr.onload = () => {
      if (xhr.status >= 200 && xhr.status < 300) resolve(xhr.getResponseHeader('ETag'))
      else reject(new Error(`Uploading failed (${xhr.status}): ${xhr.responseText.slice(0, 200)}`))
    }
    xhr.onerror = () => reject(new NetworkError('Uploading failed'))
    xhr.send(body)
  })
}

export async function uploadFile({ trpc, file, session_id, onProgress, onVerifying }: {
  trpc: UploadClient,
  file: File,
  session_id?: string,
  /** fraction of the file uploaded */
  onProgress?: (progress: number) => void,
  /** called once uploaded, while the server verifies the file */
  onVerifying?: () => void,
}): Promise<UploadFileResponse> {
  const context = { session_id }
  const prefix = session_id ? `/api/socket/${session_id}` : ''
  const { token, partSize, partCount, etagsRequired } = await trpc.components.file.uploadStart.mutate({ filename: file.name, size: file.size }, { context })
  let id: string
  try {
    const loaded = new Map<number, number>()
    const report = () => onProgress?.(file.size ? [...loaded.values()].reduce((a, b) => a + b, 0) / file.size : 1)

    // part urls are requested in batches, as they're needed (presigned urls expire)
    type PartTarget = Awaited<ReturnType<UploadClient['components']['file']['uploadPartUrls']['mutate']>>[number]
    const targets = new Map<number, Promise<PartTarget>>()
    const target = (partNumber: number) => {
      if (!targets.has(partNumber)) {
        const batch = Array.from({ length: Math.min(20, partCount - partNumber + 1) }, (_, i) => partNumber + i)
        const urls = trpc.components.file.uploadPartUrls.mutate({ token, partNumbers: batch }, { context })
        batch.forEach(n => targets.set(n, urls.then(urls => urls.find(url => url.partNumber === n)!)))
      }
      return targets.get(partNumber)!
    }
    // once a direct upload fails at the network level, the rest go through the app
    let throughApp = false
    const uploadPart = async (partNumber: number, body: ArrayBuffer) => {
      const progress = (n: number) => { loaded.set(partNumber, n); report() }
      const part = throughApp ? { partNumber, proxy: true as const } : await target(partNumber)
      if (!part.proxy) {
        try {
          const etag = await put(part.url, part.headers, body, progress)
          // an etag the browser isn't allowed to read (CORS) is fine unless the backend needs it
          if (etag || !etagsRequired) return etag ?? undefined
        } catch (e) {
          if (!(e instanceof NetworkError)) throw e
        }
        throughApp = true
      }
      for (let attempt = 1; ; attempt++) {
        try {
          return await put(`${prefix}/api/file/upload/part?token=${encodeURIComponent(token)}&part=${partNumber}`, {}, body, progress) ?? undefined
        } catch (e) {
          if (attempt >= ATTEMPTS) throw e
          await new Promise(resolve => setTimeout(resolve, 1000 * attempt))
        }
      }
    }

    const hasher = await createSHA256()
    hasher.init()
    const parts: { partNumber: number, etag?: string }[] = []
    const inflight = new Set<Promise<void>>()
    for (let partNumber = 1; partNumber <= partCount; partNumber++) {
      const start = (partNumber - 1) * partSize
      // parts are read (& hashed) in order, and uploaded a few at a time
      const body = await file.slice(start, Math.min(start + partSize, file.size)).arrayBuffer()
      hasher.update(new Uint8Array(body))
      const task = uploadPart(partNumber, body).then(etag => { parts.push({ partNumber, etag }) })
      inflight.add(task)
      task.then(() => inflight.delete(task), () => {})
      if (inflight.size >= CONCURRENCY) await Promise.race(inflight)
    }
    await Promise.all(inflight)
    report()
    ;({ id } = await trpc.components.file.uploadComplete.mutate({ token, parts, sha256: hasher.digest('hex') }, { context }))
  } catch (e) {
    await trpc.components.file.uploadAbort.mutate({ token }, { context }).catch(() => {})
    throw e
  }
  onVerifying?.()
  return await new Promise<UploadFileResponse>((resolve, reject) => {
    const subscription = trpc.components.file.uploadVerified.subscribe({ id }, {
      context,
      onData: (record) => { subscription.unsubscribe(); resolve(record) },
      onError: (error) => { subscription.unsubscribe(); reject(error) },
    })
  })
}
