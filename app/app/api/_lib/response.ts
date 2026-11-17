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

export function contentTypeOf(filename: string) {
  const ext = filename.includes('.') ? filename.slice(filename.lastIndexOf('.')).toLowerCase() : ''
  return contentTypes[ext] ?? 'application/octet-stream'
}

/**
 * The types which may be shown in the page (e.g. a PDF in an iframe) rather than downloaded:
 *  files are uploaded by users, so not anything which could run scripts (e.g. HTML, SVG).
 */
const inlineContentTypes = new Set(['application/pdf', 'image/png', 'image/jpeg', 'image/gif', 'text/plain'])

export type Disposition = 'attachment' | 'inline'

/** `Content-Disposition` for a file named `filename`, inline only for the types which may be */
export function contentDisposition(filename: string, disposition: Disposition = 'attachment') {
  const type = disposition === 'inline' && inlineContentTypes.has(contentTypeOf(filename)) ? 'inline' : 'attachment'
  // the quoted filename is a plain ascii fallback, filename* preserves the original
  return `${type}; filename="${filename.replace(/[^\x20-\x7e]|["\\]/g, '_')}"; filename*=UTF-8''${encodeURIComponent(filename)}`
}

/** Headers for downloading a file named `filename` (or showing it, see `contentDisposition`) */
export function attachmentHeaders(filename: string, size?: number, disposition: Disposition = 'attachment'): Record<string, string> {
  return {
    'Content-Type': contentTypeOf(filename),
    'Content-Disposition': contentDisposition(filename, disposition),
    'X-Content-Type-Options': 'nosniff',
    ...(size !== undefined ? { 'Content-Length': String(size) } : {}),
  }
}
