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

const contentTypes: Record<string, string> = {
  '.zip': 'application/zip',
  '.pdf': 'application/pdf',
  '.json': 'application/json',
  '.tex': 'application/x-tex',
  '.txt': 'text/plain',
  '.tsv': 'text/tab-separated-values',
  '.csv': 'text/csv',
}

/** Headers for downloading a file named `filename` */
export function attachmentHeaders(filename: string, size?: number): Record<string, string> {
  const ext = filename.includes('.') ? filename.slice(filename.lastIndexOf('.')).toLowerCase() : ''
  return {
    'Content-Type': contentTypes[ext] ?? 'application/octet-stream',
    // the quoted filename is a plain ascii fallback, filename* preserves the original
    'Content-Disposition': `attachment; filename="${filename.replace(/[^\x20-\x7e]|["\\]/g, '_')}"; filename*=UTF-8''${encodeURIComponent(filename)}`,
    ...(size !== undefined ? { 'Content-Length': String(size) } : {}),
  }
}
