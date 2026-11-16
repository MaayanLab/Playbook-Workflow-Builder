import fs from 'fs'
import { fileURLToPath } from 'url'
import { resolveDrsObject } from './drs'
import { NotFoundError, ResponseCodedError } from '@/spec/error'
import { toReadable } from '@/utils/readable'
import type { Readable } from 'stream'
import storage, { storageKey, type ByteRange } from './storage'

/** drs://{host}/{id} urls served by this instance */
function localDrsObjectId(url: string) {
  const m = /^drs:\/\/([^/]+)\/(.+)$/.exec(url)
  if (!m) return undefined
  const publicHost = process.env.PUBLIC_URL ? new URL(process.env.PUBLIC_URL).host : undefined
  return m[1] === publicHost ? m[2] : undefined
}

async function fetchStream(url: string, headers: Record<string, string> = {}, range?: ByteRange) {
  const req = await fetch(url, {
    headers: { ...headers, ...(range ? { Range: `bytes=${range.start}-${range.end ?? ''}` } : {}) },
  })
  if (req.status === 404) throw new NotFoundError()
  if (!req.ok || !req.body) throw new ResponseCodedError(req.status, await req.text())
  return toReadable(req.body)
}

/**
 * Resolve a drs url on another DRS server to its (https) access url
 */
async function resolveRemoteDrs(url: string) {
  const m = /^drs:\/\/([^/]+)\/(.+)$/.exec(url)
  if (!m) throw new Error(`Invalid drs url ${JSON.stringify(url)}`)
  const [, host, id] = m
  // like ufs: a host with an explicit port (i.e. a development server) is served over http
  const objectUrl = `${host.includes(':') ? 'http' : 'https'}://${host}/ga4gh/drs/v1/objects/${id}`
  const req = await fetch(objectUrl)
  if (req.status === 404) throw new NotFoundError()
  if (!req.ok) throw new ResponseCodedError(req.status, await req.text())
  const object = await req.json() as { access_methods?: { type: string, access_id?: string, access_url?: string | { url: string, headers?: string[] } }[] }
  const methods = (object.access_methods ?? []).filter(m => m.type === 'https' || m.type === 'http')
  const asHeaders = (headers?: string[]) => Object.fromEntries((headers ?? []).map(h => h.split(': ', 2) as [string, string]))
  for (const method of methods) {
    if (!method.access_id) continue
    const access = await fetch(`${objectUrl}/access/${method.access_id}`)
    if (access.ok) {
      const { url, headers } = await access.json() as { url: string, headers?: string[] }
      return { url, headers: asHeaders(headers) }
    }
  }
  for (const method of methods) {
    if (typeof method.access_url === 'string') return { url: method.access_url, headers: {} }
    if (method.access_url) return { url: method.access_url.url, headers: asHeaders(method.access_url.headers) }
  }
  throw new ResponseCodedError(422, `DRS object ${url} has no supported access method`)
}

/**
 * Read a file referenced by url (`file://`, `storage://`, `http(s)://` or `drs://`),
 *  optionally only the (inclusive) byte range of it
 */
export async function fileAsStream(file: { url: string }, range?: ByteRange): Promise<Readable> {
  const url = file.url
  if (url.startsWith('file://')) {
    return fs.createReadStream(fileURLToPath(url), range)
  } else if (url.startsWith('storage://')) {
    return await storage().read(storageKey(url), range)
  } else if (url.startsWith('http://') || url.startsWith('https://')) {
    return await fetchStream(url, {}, range)
  } else if (url.startsWith('drs://')) {
    const id = localDrsObjectId(url)
    if (id !== undefined) {
      const upload = await resolveDrsObject(id)
      if (!upload) throw new NotFoundError()
      return await fileAsStream(upload, range)
    }
    const access = await resolveRemoteDrs(url)
    return await fetchStream(access.url, access.headers, range)
  }
  throw new ResponseCodedError(422, `Unsupported file url ${JSON.stringify(url)}`)
}
