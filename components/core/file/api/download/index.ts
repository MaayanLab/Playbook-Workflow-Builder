import fs from 'fs'
import stream from 'stream'
import { NotFoundError, ResponseCodedError, UnauthorizedError } from '@/spec/error'
import type { Writable, Readable } from 'stream'

async function consumeReader(reader: ReadableStreamDefaultReader, writeable: Writable) {
  const { value, done } = await reader.read()
  if (value) writeable.write(value)
  if (!done) process.nextTick(consumeReader, reader, writeable)
}

function readerToReadable(reader: ReadableStream): Readable {
  const duplex = new stream.Duplex()
  consumeReader(reader.getReader(), duplex)
  return duplex
}

export async function downloadFile(file: { url: string }): Promise<Readable> {
  if (file.url.startsWith('file://')) {
    if (process.env.NODE_ENV !== 'development') throw new UnauthorizedError()
    return fs.createReadStream(file.url.slice('file://'.length))
  } else {
    const req = await fetch(file.url)
    if (req.status === 404) throw new NotFoundError()
    else if (req.body) {
      return readerToReadable(req.body)
    }
    throw new ResponseCodedError(req.status, await req.text())
  }
}
