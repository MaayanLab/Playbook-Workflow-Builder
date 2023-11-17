import { Readable } from 'stream'

async function *yieldReadableStream(stream: ReadableStream<Uint8Array>) {
  const reader = stream.getReader()
  while (true) {
    const data = await reader.read()
    if (data.value) yield data.value
    if (data.done) break
  }
}

export function toReadable(stream: Readable | ReadableStream<Uint8Array>) {
  if ('pipe' in stream) return stream
  return Readable.from(yieldReadableStream(stream))
}
