import fs from 'fs'
import path from 'path'
import { NotFoundError, ResponseCodedError, UnauthorizedError } from '@/spec/error'
import { toReadable } from '@/utils/readable'
import type { Readable } from 'stream'

export async function fileAsStream(file: { url: string }): Promise<Readable> {
  let url = file.url
  if (url.startsWith((process.env.PUBLIC_URL||'').replace('https?:/', 'drs:/'))) {
    url = url.replace(
      /^drs:\/\/([^\/]+)\/([^\/]+)$/, 
      'https://$1/ga4gh/drs/v1/objects/$2/access/https'
    )
  }
  if (url.startsWith('file://')) {
    // TODO: this should be enabled when we have DRS support in prod
    // if (process.env.NODE_ENV !== 'development') throw new UnauthorizedError()
    return fs.createReadStream(url.slice('file://'.length).split(path.posix.sep).join(path.sep))
  } else {
    const req = await fetch(url)
    if (req.status === 404) throw new NotFoundError()
    else if (req.body) {
      return toReadable(req.body)
    }
    throw new ResponseCodedError(req.status, await req.text())
  }
}
