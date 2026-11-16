import db from '@/app/db'
import { getRequestSessionWithId } from '@/app/extensions/auth/server'
import { forward } from '@/app/extensions/cavatica/tunnel'

const hopByHop = new Set(['connection', 'keep-alive', 'transfer-encoding', 'upgrade', 'host', 'content-length', 'content-encoding', 'accept-encoding'])
const HEAD_TIMEOUT_MS = 30_000

async function handler(request: Request, props: { params: Promise<{ id: string, path?: string[] }> }) {
  const params = await props.params
  const session = await getRequestSessionWithId(request)
  if (!session?.user) return new Response('Unauthorized', { status: 401 })
  const proxy = await db.selectFrom('proxy_session')
    .innerJoin('user_integrations', 'user_integrations.id', 'proxy_session.user')
    .select('user_integrations.cavatica_api_key')
    .where('proxy_session.id', '=', params.id)
    .where('proxy_session.user', '=', session.user.id)
    .executeTakeFirst()
  if (!proxy?.cavatica_api_key) return new Response('Session not found', { status: 404 })

  const upstreamPath = `/${(params.path ?? []).map(encodeURIComponent).join('/')}${new URL(request.url).search}`
  const body = ['GET', 'HEAD'].includes(request.method) ? undefined : Buffer.from(await request.arrayBuffer()).toString('base64')
  let status = 200
  let responseHeaders: Record<string, string> = {}
  let start!: () => void
  const started = new Promise<void>(resolve => { start = resolve })
  let abort = () => {}
  let timer: ReturnType<typeof setTimeout>

  const stream = new ReadableStream<Uint8Array>({
    start(controller) {
      timer = setTimeout(() => {
        status = 504
        abort()
        controller.enqueue(new TextEncoder().encode('Worker did not respond'))
        controller.close()
        start()
      }, HEAD_TIMEOUT_MS)
      try {
        abort = forward(params.id, {
          method: request.method,
          path: upstreamPath,
          headers: {
            ...Object.fromEntries(request.headers.entries()),
            authorization: `Token ${proxy.cavatica_api_key}`,
          },
          body,
        }, {
          onHead(code, headers) {
            clearTimeout(timer)
            status = code
            responseHeaders = Object.fromEntries(Object.entries(headers).filter(([key]) => !hopByHop.has(key.toLowerCase())))
            start()
          },
          onChunk(data) { controller.enqueue(typeof data === 'string' ? new TextEncoder().encode(data) : data) },
          onEnd() { clearTimeout(timer); controller.close(); start() },
          onError(error) {
            clearTimeout(timer)
            status = 502
            controller.enqueue(new TextEncoder().encode(error.message))
            controller.close()
            start()
          },
        })
      } catch {
        clearTimeout(timer)
        status = 503
        controller.enqueue(new TextEncoder().encode('Worker not connected'))
        controller.close()
        start()
      }
      request.signal.addEventListener('abort', abort, { once: true })
    },
    cancel() { clearTimeout(timer); abort() },
  })
  await started
  return new Response(request.method === 'HEAD' ? null : stream, { status, headers: responseHeaders })
}

export { handler as GET, handler as POST, handler as PUT, handler as PATCH, handler as DELETE, handler as HEAD }
