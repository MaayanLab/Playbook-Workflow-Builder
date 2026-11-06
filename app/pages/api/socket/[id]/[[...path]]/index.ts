/**
 * Forward a request to the CAVATICA worker of a session through the tunnel
 *  (see @/app/extensions/cavatica/tunnel), streaming the response back.
 */
import type { NextApiRequest, NextApiResponse } from 'next'
import getRawBody from 'raw-body'
import db from '@/app/db'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { forward } from '@/app/extensions/cavatica/tunnel'
import * as dict from '@/utils/dict'

export const config = {
  api: {
    bodyParser: false,
    // responses are streamed (and may be long lived, e.g. subscriptions)
    responseLimit: false,
  },
}

/** headers which only apply to a single connection */
const hopByHop = new Set(['connection', 'keep-alive', 'transfer-encoding', 'upgrade', 'host', 'content-length', 'content-encoding', 'accept-encoding'])
const forwardable = (headers: Record<string, string>) => dict.filter(headers, ({ key }) => !hopByHop.has(key.toLowerCase()))

/** time to wait for the worker to start responding */
const HEAD_TIMEOUT_MS = 30_000

export default async function Forward(req: NextApiRequest, res: NextApiResponse) {
  const session = await getServerSessionWithId(req, res)
  if (!session?.user) return void res.status(401).end('Unauthorized')
  const session_id = req.query.id as string
  const proxy_session = await db.selectFrom('proxy_session')
    .innerJoin('user_integrations', 'user_integrations.id', 'proxy_session.user')
    .select('user_integrations.cavatica_api_key')
    .where('proxy_session.id', '=', session_id)
    .where('proxy_session.user', '=', session.user.id)
    .executeTakeFirst()
  if (!proxy_session?.cavatica_api_key) return void res.status(404).end('Session not found')
  const path = `/${((req.query.path ?? []) as string[]).map(encodeURIComponent).join('/')}`
  const search = req.url?.includes('?') ? req.url.slice(req.url.indexOf('?')) : ''
  const body = ['HEAD', 'GET'].includes(req.method ?? '') ? undefined : (await getRawBody(req)).toString('base64')
  await new Promise<void>((resolve) => {
    let started = false
    const headTimeout = setTimeout(() => {
      abort()
      if (!res.headersSent) res.status(504).end('Worker did not respond')
      resolve()
    }, HEAD_TIMEOUT_MS)
    let abort = () => {}
    try {
      abort = forward(session_id, {
        method: req.method ?? 'GET',
        path: `${path}${search}`,
        headers: {
          // only accept, content-type & authorization are forwarded
          ...dict.init(dict.items(dict.fromIncomingHeaders(req.headers)).flatMap(({ key, value }) =>
            value === undefined ? [] : [{ key: key.toString(), value: Array.isArray(value) ? value.join(', ') : value }])),
          // the worker treats requests with its own key as its user's
          authorization: `Token ${proxy_session.cavatica_api_key}`,
        },
        body,
      }, {
        onHead: (status, headers) => {
          clearTimeout(headTimeout)
          started = true
          res.writeHead(status, forwardable(headers))
          res.flushHeaders()
        },
        onChunk: (data) => { res.write(data) },
        onEnd: () => { res.end(); resolve() },
        onError: (error) => {
          clearTimeout(headTimeout)
          if (!started) res.status(502).end(error.message)
          else res.end()
          resolve()
        },
      })
    } catch (error) {
      clearTimeout(headTimeout)
      res.status(503).end('Worker not connected')
      return resolve()
    }
    // the client went away, stop the request on the worker
    res.on('close', () => {
      clearTimeout(headTimeout)
      abort()
      resolve()
    })
  })
}
