import fs from 'fs'
import path from 'path'
import { NotFoundError, ResponseCodedError, UnauthorizedError } from '@/spec/error'
import { toReadable } from '@/utils/readable'
import type { Readable } from 'stream'
import { pythonStream } from '@/utils/python'

export async function fileAsStream(file: { url: string }): Promise<Readable> {
  let url = file.url
  if (url.startsWith('file://')) {
    return fs.createReadStream(url.slice('file://'.length).split(path.posix.sep).join(path.sep))
  } else if (url.startsWith('http://') || url.startsWith('https://')) {
    const req = await fetch(url)
    if (req.status === 404) throw new NotFoundError()
    else if (req.body) {
      return toReadable(req.body)
    }
    throw new ResponseCodedError(req.status, await req.text())
  } else {
    return pythonStream('components.core.file.file_read_stream', {
      kargs: [file],
    })
  }
}
