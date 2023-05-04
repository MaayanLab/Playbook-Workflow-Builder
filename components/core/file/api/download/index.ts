import fs from 'fs'
import { NotFoundError, ResponseCodedError, UnauthorizedError } from '@/spec/error'
import { toReadable } from '@/utils/readable'
import type { Readable } from 'stream'

export async function fileAsStream(file: { url: string }): Promise<Readable> {
  if (file.url.startsWith('file://')) {
    if (process.env.NODE_ENV !== 'development') throw new UnauthorizedError()
    return fs.createReadStream(file.url.slice('file://'.length))
  } else {
    const req = await fetch(file.url)
    if (req.status === 404) throw new NotFoundError()
    else if (req.body) {
      return toReadable(req.body)
    }
    throw new ResponseCodedError(req.status, await req.text())
  }
}
