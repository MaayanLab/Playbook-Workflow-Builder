/**
 * Content addressed file storage: files are stored by their sha256 (`storage://<sha256>`).
 *
 * The backend is selected by the environment:
 * - S3 (`S3_BUCKET`, optionally `S3_PREFIX`, `S3_ENDPOINT` for S3 compatible services & `S3_REGION`,
 *    credentials come from the standard AWS environment (`AWS_ACCESS_KEY_ID`, ...)).
 *    A `UFS_STORAGE` configuration of an S3 bucket is also understood, and maps to the same
 *    objects ufs stored, so existing files remain accessible.
 * - A CAVATICA project (`CAVATICA_STORAGE_PROJECT`, used by CAVATICA workers) through the
 *    Seven Bridges API client.
 * - The local filesystem otherwise (`STORAGE_PATH`, default `data/ufs` where ufs stored them).
 */
import fs from 'fs'
import path from 'path'
import { Readable } from 'stream'
import { pipeline } from 'stream/promises'
import {
  S3Client, GetObjectCommand, HeadObjectCommand, DeleteObjectCommand,
  CreateMultipartUploadCommand, UploadPartCommand, CompleteMultipartUploadCommand, AbortMultipartUploadCommand, ListPartsCommand,
  CopyObjectCommand, UploadPartCopyCommand,
} from '@aws-sdk/client-s3'
import { Upload } from '@aws-sdk/lib-storage'
import { getSignedUrl } from '@aws-sdk/s3-request-presigner'
import python from '@/utils/python'
import { NotFoundError, ResponseCodedError } from '@/spec/error'
import cache from '@/utils/global_cache'

export type ByteRange = { start: number, end?: number }

export interface StorageBackend {
  readonly type: 's3' | 'filesystem' | 'cavatica'
  /** The size of the stored object, undefined if it doesn't exist */
  stat(key: string): Promise<{ size: number } | undefined>
  /** Read the object (or the inclusive byte range of it) */
  read(key: string, range?: ByteRange): Promise<Readable>
  /** Move a local file into storage (it's removed from its local path) */
  putFile(key: string, localPath: string): Promise<void>
  /** A url which can be used (for a limited time) to GET the object directly (supporting Range) */
  presignedGet?(key: string, opts: { filename?: string, contentType?: string, expiresIn?: number }): Promise<string>
  delete(key: string): Promise<void>
  /** Copy an object within storage (without transferring it through the app) */
  copy(srcKey: string, dstKey: string): Promise<void>

  /**
   * Multipart uploads, used to upload files from the browser in parts (see ./direct-upload.ts).
   *  Parts are PUT directly to storage where the backend can give out urls for them (S3, CAVATICA),
   *  or sent through the app otherwise (or when PUTting them directly fails).
   */
  createMultipart(key: string, size: number): Promise<{ uploadId: string, partSize: number }>
  /** A url (& headers) to PUT a part to directly, undefined if parts must go through `uploadPart` */
  partUrl(key: string, uploadId: string, partNumber: number): Promise<{ url: string, headers: Record<string, string> } | undefined>
  /** Upload a part through the app, returning its etag */
  uploadPart(key: string, uploadId: string, partNumber: number, body: Readable, size: number): Promise<string>
  /** Assemble the object from its parts (a part's etag may be missing where the backend can look it up) */
  completeMultipart(key: string, uploadId: string, parts: MultipartPart[]): Promise<void>
  abortMultipart(key: string, uploadId: string): Promise<void>
}

export type MultipartPart = { partNumber: number, etag?: string }

const MiB = 1024 * 1024
/** S3 (& S3 based services) limit multipart uploads to 10000 parts of at least 5MiB */
export const partSizeFor = (size: number, minimum: number) => Math.max(minimum, Math.ceil(size / 10000 / MiB) * MiB)

/** PUT a part to a (presigned) url, returning its etag */
async function putPart(url: string, headers: Record<string, string>, body: Readable, size: number) {
  const res = await fetch(url, {
    method: 'PUT',
    headers: { ...headers, 'Content-Length': String(size) },
    body: Readable.toWeb(body) as ReadableStream,
    // @ts-ignore -- required by node to stream a request body
    duplex: 'half',
  })
  if (!res.ok) throw new ResponseCodedError(res.status, await res.text())
  const etag = res.headers.get('etag')
  if (!etag) throw new Error('Part upload returned no etag')
  return etag
}

const uploadIdC = /^[A-Za-z0-9._~-]+$/
function checkUploadId(uploadId: string) {
  if (!uploadIdC.test(uploadId)) throw new Error(`Invalid upload id ${JSON.stringify(uploadId)}`)
  return uploadId
}

/** `storage://<key>` => `<key>` */
export function storageKey(url: string) {
  const m = /^storage:\/\/([A-Za-z0-9._-]+)$/.exec(url)
  if (!m) throw new Error(`Invalid storage url ${JSON.stringify(url)}`)
  return m[1]
}

async function moveFile(src: string, dst: string) {
  await fs.promises.mkdir(path.dirname(dst), { recursive: true })
  try {
    await fs.promises.rename(src, dst)
  } catch (e: any) {
    // across devices (e.g. from the os tmpdir), copy instead
    if (e?.code !== 'EXDEV') throw e
    await fs.promises.copyFile(src, dst)
    await fs.promises.unlink(src)
  }
}

export class FilesystemStorage implements StorageBackend {
  readonly type = 'filesystem'
  constructor(private root: string) {}
  private path = (key: string) => path.join(this.root, storageKey(`storage://${key}`))
  stat = async (key: string) => {
    const stat = await fs.promises.stat(this.path(key)).catch(() => undefined)
    return stat?.isFile() ? { size: stat.size } : undefined
  }
  read = async (key: string, range?: ByteRange) => fs.createReadStream(this.path(key), range)
  putFile = async (key: string, localPath: string) => {
    // content addressed: if we already have it, we don't need another copy
    if (await this.stat(key)) await fs.promises.unlink(localPath)
    else await moveFile(localPath, this.path(key))
  }
  delete = async (key: string) => { await fs.promises.rm(this.path(key), { force: true }) }
  copy = async (srcKey: string, dstKey: string) => {
    // a hard link doesn't take up any more space
    await fs.promises.link(this.path(srcKey), this.path(dstKey)).catch(async (e: any) => {
      if (e?.code === 'EEXIST') return
      await fs.promises.copyFile(this.path(srcKey), this.path(dstKey))
    })
  }

  // parts are staged in a directory, then concatenated
  private staging = (uploadId: string) => path.join(this.root, '.multipart', checkUploadId(uploadId))
  createMultipart = async (_key: string, size: number) => {
    await this.sweepAbandoned()
    const uploadId = crypto.randomUUID()
    await fs.promises.mkdir(this.staging(uploadId), { recursive: true })
    // parts go through the app, so they're kept small enough to upload quickly
    return { uploadId, partSize: partSizeFor(size, 8 * MiB) }
  }
  /** remove the parts of uploads which were never completed (upload tokens are valid for a day) */
  private sweepAbandoned = async () => {
    const root = path.join(this.root, '.multipart')
    const entries = await fs.promises.readdir(root).catch(() => [] as string[])
    for (const entry of entries) {
      const stat = await fs.promises.stat(path.join(root, entry)).catch(() => undefined)
      if (stat && Date.now() - stat.mtimeMs > 2 * 24 * 60 * 60 * 1000) {
        await fs.promises.rm(path.join(root, entry), { recursive: true, force: true })
      }
    }
  }
  partUrl = async () => undefined
  uploadPart = async (_key: string, uploadId: string, partNumber: number, body: Readable, size: number) => {
    const partPath = path.join(this.staging(uploadId), String(partNumber))
    await pipeline(body, fs.createWriteStream(partPath))
    const { size: written } = await fs.promises.stat(partPath)
    if (written !== size) throw new ResponseCodedError(400, `Part ${partNumber} is ${written} bytes, expected ${size}`)
    return `"${partNumber}-${size}"`
  }
  completeMultipart = async (key: string, uploadId: string, parts: MultipartPart[]) => {
    const staging = this.staging(uploadId)
    const assembled = path.join(staging, 'assembled')
    const out = fs.createWriteStream(assembled)
    for (const { partNumber } of [...parts].sort((a, b) => a.partNumber - b.partNumber)) {
      await pipeline(fs.createReadStream(path.join(staging, String(partNumber))), out, { end: false })
    }
    await new Promise<void>((resolve, reject) => out.end((err?: Error | null) => err ? reject(err) : resolve()))
    await moveFile(assembled, this.path(key))
    await fs.promises.rm(staging, { recursive: true, force: true })
  }
  abortMultipart = async (_key: string, uploadId: string) => {
    await fs.promises.rm(this.staging(uploadId), { recursive: true, force: true })
  }
}

export type S3StorageConfig = {
  bucket: string,
  /** key prefix, e.g. `uploads/` */
  prefix: string,
  region: string,
  endpoint?: string,
  credentials?: { accessKeyId: string, secretAccessKey: string },
}

export class S3Storage implements StorageBackend {
  readonly type = 's3'
  private client: S3Client
  constructor(private config: S3StorageConfig) {
    this.client = new S3Client({
      region: config.region,
      endpoint: config.endpoint,
      // S3 compatible services (minio, ...) are typically addressed by path
      forcePathStyle: config.endpoint !== undefined && !/(^|\.)amazonaws\.com(:\d+)?$/.test(new URL(config.endpoint).host),
      credentials: config.credentials,
      // by default the SDK adds checksums to requests, presigned upload urls would then carry
      //  the checksum of an empty body (& S3 compatible services don't all support them)
      requestChecksumCalculation: 'WHEN_REQUIRED',
      responseChecksumValidation: 'WHEN_REQUIRED',
    })
  }
  private key = (key: string) => `${this.config.prefix}${storageKey(`storage://${key}`)}`
  stat = async (key: string) => {
    try {
      const head = await this.client.send(new HeadObjectCommand({ Bucket: this.config.bucket, Key: this.key(key) }))
      return { size: head.ContentLength ?? 0 }
    } catch (e: any) {
      if (e?.$metadata?.httpStatusCode === 404 || e?.name === 'NotFound') return undefined
      throw e
    }
  }
  read = async (key: string, range?: ByteRange) => {
    const res = await this.client.send(new GetObjectCommand({
      Bucket: this.config.bucket,
      Key: this.key(key),
      Range: range ? `bytes=${range.start}-${range.end ?? ''}` : undefined,
    }))
    return res.Body as Readable
  }
  putFile = async (key: string, localPath: string) => {
    if (!(await this.stat(key))) {
      // multipart for large files
      await new Upload({
        client: this.client,
        params: { Bucket: this.config.bucket, Key: this.key(key), Body: fs.createReadStream(localPath) },
      }).done()
    }
    await fs.promises.unlink(localPath)
  }
  presignedGet = async (key: string, { filename, contentType, expiresIn = 60 * 60 }: { filename?: string, contentType?: string, expiresIn?: number }) => {
    return await getSignedUrl(this.client, new GetObjectCommand({
      Bucket: this.config.bucket,
      Key: this.key(key),
      ResponseContentDisposition: filename ? `attachment; filename="${filename.replace(/[^\x20-\x7e]|["\\]/g, '_')}"; filename*=UTF-8''${encodeURIComponent(filename)}` : undefined,
      ResponseContentType: contentType,
    }), { expiresIn })
  }
  delete = async (key: string) => {
    await this.client.send(new DeleteObjectCommand({ Bucket: this.config.bucket, Key: this.key(key) }))
  }
  copy = async (srcKey: string, dstKey: string) => {
    const CopySource = `${this.config.bucket}/${this.key(srcKey)}`
    const size = (await this.stat(srcKey))?.size
    if (size === undefined) throw new NotFoundError()
    if (size <= 5 * 1024 * MiB) {
      await this.client.send(new CopyObjectCommand({ Bucket: this.config.bucket, Key: this.key(dstKey), CopySource }))
      return
    }
    // larger objects are copied in parts
    const { uploadId, partSize } = await this.createMultipart(dstKey, size)
    try {
      const parts: MultipartPart[] = []
      for (let partNumber = 1, start = 0; start < size; partNumber++, start += partSize) {
        const res = await this.client.send(new UploadPartCopyCommand({
          Bucket: this.config.bucket, Key: this.key(dstKey), UploadId: uploadId, PartNumber: partNumber,
          CopySource, CopySourceRange: `bytes=${start}-${Math.min(start + partSize, size) - 1}`,
        }))
        parts.push({ partNumber, etag: res.CopyPartResult?.ETag })
      }
      await this.completeMultipart(dstKey, uploadId, parts)
    } catch (e) {
      await this.abortMultipart(dstKey, uploadId).catch(() => {})
      throw e
    }
  }

  createMultipart = async (key: string, size: number) => {
    const res = await this.client.send(new CreateMultipartUploadCommand({ Bucket: this.config.bucket, Key: this.key(key) }))
    if (!res.UploadId) throw new Error('S3 returned no upload id')
    return { uploadId: res.UploadId, partSize: partSizeFor(size, 16 * MiB) }
  }
  partUrl = async (key: string, uploadId: string, partNumber: number) => ({
    url: await getSignedUrl(this.client, new UploadPartCommand({
      Bucket: this.config.bucket, Key: this.key(key), UploadId: uploadId, PartNumber: partNumber,
    }), { expiresIn: 60 * 60 }),
    headers: {},
  })
  uploadPart = async (key: string, uploadId: string, partNumber: number, body: Readable, size: number) => {
    const res = await this.client.send(new UploadPartCommand({
      Bucket: this.config.bucket, Key: this.key(key), UploadId: uploadId, PartNumber: partNumber,
      Body: body, ContentLength: size,
    }))
    if (!res.ETag) throw new Error('S3 returned no etag')
    return res.ETag
  }
  completeMultipart = async (key: string, uploadId: string, parts: MultipartPart[]) => {
    // the browser can only read etags if the bucket's CORS configuration exposes them, so
    //  the uploaded parts are looked up instead
    const etags = new Map<number, string>()
    let marker: string | undefined
    do {
      const res = await this.client.send(new ListPartsCommand({ Bucket: this.config.bucket, Key: this.key(key), UploadId: uploadId, PartNumberMarker: marker }))
      for (const part of res.Parts ?? []) if (part.PartNumber && part.ETag) etags.set(part.PartNumber, part.ETag)
      marker = res.IsTruncated ? res.NextPartNumberMarker : undefined
    } while (marker)
    await this.client.send(new CompleteMultipartUploadCommand({
      Bucket: this.config.bucket, Key: this.key(key), UploadId: uploadId,
      MultipartUpload: {
        Parts: [...parts].sort((a, b) => a.partNumber - b.partNumber).map(({ partNumber, etag }) => {
          const ETag = etags.get(partNumber) ?? etag
          if (!ETag) throw new ResponseCodedError(400, `Part ${partNumber} wasn't uploaded`)
          return { PartNumber: partNumber, ETag }
        }),
      },
    }))
  }
  abortMultipart = async (key: string, uploadId: string) => {
    await this.client.send(new AbortMultipartUploadCommand({ Bucket: this.config.bucket, Key: this.key(key), UploadId: uploadId }))
  }
}

/**
 * A CAVATICA project (`CAVATICA_STORAGE_PROJECT`), through the Seven Bridges API client in python
 *  (see components/core/file/cavatica_storage.py). Objects are read from their (presigned) download urls.
 */
export class CavaticaStorage implements StorageBackend {
  readonly type = 'cavatica'
  stat = async (key: string) => await python<{ size: number } | null>('components.core.file.cavatica_storage.stat', { kargs: [storageKey(`storage://${key}`)] }) ?? undefined
  // download urls are reused for a few minutes (well within their expiry), reading an object in
  //  chunks would otherwise ask the API for one per chunk (it's limited to 1000 requests / 5 min)
  private downloadUrls = new Map<string, { url: string, expires: number }>()
  private downloadUrl = async (key: string) => {
    const cached = this.downloadUrls.get(key)
    if (cached && cached.expires > Date.now()) return cached.url
    const url = await python<string | null>('components.core.file.cavatica_storage.download_url', { kargs: [storageKey(`storage://${key}`)] })
    if (!url) throw new NotFoundError()
    this.downloadUrls.set(key, { url, expires: Date.now() + 5 * 60 * 1000 })
    return url
  }
  read = async (key: string, range?: ByteRange) => {
    const res = await fetch(await this.downloadUrl(key), {
      headers: range ? { Range: `bytes=${range.start}-${range.end ?? ''}` } : {},
    })
    if (!res.ok || !res.body) throw new ResponseCodedError(res.status, await res.text())
    return Readable.fromWeb(res.body as any)
  }
  putFile = async (key: string, localPath: string) => {
    await python('components.core.file.cavatica_storage.put', { kargs: [storageKey(`storage://${key}`), localPath] })
  }
  // (the platform controls the downloaded file's name)
  presignedGet = async (key: string) => await this.downloadUrl(key)
  delete = async (key: string) => { await python('components.core.file.cavatica_storage.delete', { kargs: [storageKey(`storage://${key}`)] }) }
  copy = async (srcKey: string, dstKey: string) => {
    await python('components.core.file.cavatica_storage.copy', { kargs: [storageKey(`storage://${srcKey}`), storageKey(`storage://${dstKey}`)] })
  }

  createMultipart = async (key: string, size: number) => {
    const { upload_id, part_size } = await python<{ upload_id: string, part_size: number }>('components.core.file.cavatica_storage.create_multipart', {
      kargs: [storageKey(`storage://${key}`), size, partSizeFor(size, 16 * MiB)],
    })
    return { uploadId: upload_id, partSize: part_size }
  }
  partUrl = async (_key: string, uploadId: string, partNumber: number) =>
    await python<{ url: string, headers: Record<string, string> }>('components.core.file.cavatica_storage.part_url', { kargs: [checkUploadId(uploadId), partNumber] })
  uploadPart = async (key: string, uploadId: string, partNumber: number, body: Readable, size: number) => {
    const { url, headers } = await this.partUrl(key, uploadId, partNumber)
    return await putPart(url, headers, body, size)
  }
  completeMultipart = async (_key: string, uploadId: string, parts: MultipartPart[]) => {
    // the platform needs every part's etag reported
    const missing = parts.find(part => !part.etag)
    if (missing) throw new ResponseCodedError(400, `Part ${missing.partNumber} has no etag`)
    await python('components.core.file.cavatica_storage.complete_multipart', {
      kargs: [checkUploadId(uploadId), parts.map(({ partNumber, etag }) => ({ part_number: partNumber, etag }))],
    })
  }
  abortMultipart = async (_key: string, uploadId: string) => {
    await python('components.core.file.cavatica_storage.abort_multipart', { kargs: [checkUploadId(uploadId)] })
  }
}

type UfsConfig = { cls: string, [key: string]: any }

/**
 * Translate a ufs configuration of an S3 bucket (possibly wrapped in a Prefix or caches)
 *  to the S3 bucket & key prefix ufs used, undefined if it isn't one.
 */
export function s3ConfigFromUfs(config: UfsConfig): S3StorageConfig | undefined {
  // ufs paths are posix paths, joined from the innermost prefix outwards
  function resolve(config: UfsConfig): { s3: UfsConfig, prefix: string } | undefined {
    if (config.cls === 'ufs.impl.s3.S3') return { s3: config, prefix: '/' }
    if (!config.ufs) return undefined
    const inner = resolve(config.ufs)
    if (!inner) return undefined
    if (config.cls === 'ufs.impl.prefix.Prefix') {
      return { s3: inner.s3, prefix: path.posix.join(inner.prefix, String(config.prefix ?? '/'), '/') }
    }
    // caches & synchronous wrappers don't change paths
    if (['ufs.impl.sync.Sync', 'ufs.impl.writecache.Writecache', 'ufs.impl.dircache.DirCache', 'ufs.impl.asyn.Async'].includes(config.cls)) return inner
    return undefined
  }
  const resolved = resolve(config)
  if (!resolved) return undefined
  // the first path component is the bucket
  const [bucket, ...rest] = resolved.prefix.split('/').filter(Boolean)
  if (!bucket) return undefined
  const { access_key, secret_access_key, endpoint_url } = resolved.s3
  return {
    bucket,
    prefix: rest.length ? `${rest.join('/')}/` : '',
    region: process.env.S3_REGION ?? process.env.AWS_REGION ?? 'us-east-1',
    endpoint: endpoint_url && endpoint_url !== 'https://s3.amazonaws.com' ? endpoint_url : undefined,
    credentials: access_key && secret_access_key ? { accessKeyId: access_key, secretAccessKey: secret_access_key } : undefined,
  }
}

export function storageFromEnv(env: NodeJS.ProcessEnv = process.env): StorageBackend {
  // (CAVATICA workers store files in the user's project, whatever else the environment configures)
  if (env.CAVATICA_STORAGE_PROJECT) return new CavaticaStorage()
  if (env.S3_BUCKET) {
    return new S3Storage({
      bucket: env.S3_BUCKET,
      prefix: env.S3_PREFIX ? `${env.S3_PREFIX.replace(/^\/+|\/+$/g, '')}/` : '',
      region: env.S3_REGION ?? env.AWS_REGION ?? 'us-east-1',
      endpoint: env.S3_ENDPOINT || undefined,
    })
  }
  if (env.UFS_STORAGE) {
    const s3 = s3ConfigFromUfs(JSON.parse(env.UFS_STORAGE))
    if (!s3) throw new Error('UFS_STORAGE is only supported for S3 buckets, configure S3_BUCKET (or CAVATICA_STORAGE_PROJECT) instead')
    return new S3Storage(s3)
  }
  return new FilesystemStorage(env.STORAGE_PATH ?? path.join(env.PYTHON_ROOT ?? env.APP_ROOT ?? process.cwd(), 'data', 'ufs'))
}

export default function storage() {
  return cache('storage', () => storageFromEnv())
}
