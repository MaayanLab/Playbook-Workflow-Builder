/**
 * Download a DRS object (an upload). Objects in storage which can be read directly (S3)
 *  are redirected to a presigned url, otherwise they're streamed (supporting Range requests).
 */
import { resolveDrsObject } from '@/components/core/file/server/drs'
import { fileAsStream } from '@/components/core/file/server/download'
import storage, { storageKey } from '@/components/core/file/server/storage'
import { NotFoundError } from '@/spec/error'
import { attachmentHeaders, nodeStreamResponse, routeResponse } from '@/app/app/api/_lib/response'

type Params = { params: Promise<{ object_id: string, access_id: string }> }

async function getUpload(props: Params) {
  const params = await props.params
  if (params.access_id !== 'https') throw new NotFoundError()
  const upload = await resolveDrsObject(params.object_id)
  if (!upload) throw new NotFoundError()
  return upload
}

/** a single `bytes=start-end` range (suffix & multiple ranges aren't supported) */
function parseRange(header: string | null, size: number) {
  const m = header ? /^bytes=(\d+)-(\d*)$/.exec(header.trim()) : null
  if (!m) return undefined
  const start = +m[1]
  const end = Math.min(m[2] ? +m[2] : size - 1, size - 1)
  if (start > end) return null
  return { start, end }
}

export async function HEAD(_request: Request, props: Params) {
  return routeResponse(async () => {
    const upload = await getUpload(props)
    return new Response(null, { headers: { ...attachmentHeaders(upload.filename, upload.size), 'Accept-Ranges': 'bytes' } })
  })
}

async function handler(request: Request, props: Params) {
  return routeResponse(async () => {
    if (request.method === 'POST') await request.json()
    const upload = await getUpload(props)
    const backend = storage()
    if (upload.url.startsWith('storage://') && backend.presignedGet) {
      // the client (incl. its Range header) is redirected to the object itself
      return Response.redirect(await backend.presignedGet(storageKey(upload.url), { filename: upload.filename }), 302)
    }
    const headers = { ...attachmentHeaders(upload.filename, upload.size), 'Accept-Ranges': 'bytes' }
    const range = parseRange(request.headers.get('range'), upload.size)
    if (range === null) {
      return new Response(null, { status: 416, headers: { 'Content-Range': `bytes */${upload.size}` } })
    } else if (range) {
      return nodeStreamResponse(await fileAsStream(upload, range), {
        status: 206,
        headers: {
          ...headers,
          'Content-Length': String(range.end - range.start + 1),
          'Content-Range': `bytes ${range.start}-${range.end}/${upload.size}`,
        },
      })
    }
    return nodeStreamResponse(await fileAsStream(upload), { headers })
  })
}

export { handler as GET, handler as POST }
