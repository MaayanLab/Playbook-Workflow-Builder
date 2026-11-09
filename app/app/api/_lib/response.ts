import { ResponseCodedError } from '@/spec/error'
import { Readable } from 'stream'

export async function routeResponse(fn: () => Promise<Response>) {
  try {
    return await fn()
  } catch (error) {
    const status = ResponseCodedError.isinstance(error) ? (error as ResponseCodedError).error_code : 500
    if (status === 500) console.error(error)
    return new Response((error as Error).message || String(error), { status })
  }
}

export function nodeStreamResponse(stream: Readable, init?: ResponseInit) {
  return new Response(Readable.toWeb(stream) as ReadableStream, init)
}
