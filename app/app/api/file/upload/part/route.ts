/**
 * A part of a direct upload sent through the app, for storage without urls to PUT parts to
 *  (the local filesystem), or when PUTting a part to storage directly failed (e.g. CORS).
 *  The upload token (see components/core/file/server/direct-upload.ts) authorizes it.
 */
import { Readable } from 'stream'
import { uploadPartThroughApp } from '@/components/core/file/server/direct-upload'
import { ResponseCodedError } from '@/spec/error'
import { routeResponse } from '@/app/app/api/_lib/response'

export const runtime = 'nodejs'

export async function PUT(request: Request) {
  return routeResponse(async () => {
    const query = new URL(request.url).searchParams
    const token = query.get('token')
    const partNumber = Number(query.get('part'))
    if (!token || !request.body) throw new ResponseCodedError(400, 'Expected a token, part number & body')
    const contentLength = request.headers.get('content-length')
    const etag = await uploadPartThroughApp(token, partNumber, Readable.fromWeb(request.body as any), contentLength ? Number(contentLength) : undefined)
    return new Response(null, { status: 200, headers: { ETag: etag, 'Access-Control-Expose-Headers': 'ETag' } })
  })
}
