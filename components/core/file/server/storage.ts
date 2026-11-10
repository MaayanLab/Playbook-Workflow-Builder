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
import { S3Client, GetObjectCommand, HeadObjectCommand, DeleteObjectCommand } from '@aws-sdk/client-s3'
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
}

/**
 * A CAVATICA project (`CAVATICA_STORAGE_PROJECT`), through the Seven Bridges API client in python
 *  (see components/core/file/cavatica_storage.py). Objects are read from their (presigned) download urls.
 */
export class CavaticaStorage implements StorageBackend {
  readonly type = 'cavatica'
  stat = async (key: string) => await python<{ size: number } | null>('components.core.file.cavatica_storage.stat', { kargs: [storageKey(`storage://${key}`)] }) ?? undefined
  private downloadUrl = async (key: string) => {
    const url = await python<string | null>('components.core.file.cavatica_storage.download_url', { kargs: [storageKey(`storage://${key}`)] })
    if (!url) throw new NotFoundError()
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
