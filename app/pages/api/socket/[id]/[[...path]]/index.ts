import type { Server } from "socket.io"
import type { NextApiRequest, NextApiResponse } from "next"
import { randomUUID } from "crypto"
import getRawBody from 'raw-body'
import { TimeoutError, UnauthorizedError } from "@/spec/error"
import { emitter } from '@/app/extensions/socket/cavatica'
import * as dict from '@/utils/dict'
import { getServerSessionWithId } from "@/app/extensions/next-auth/helpers"
import db from '@/app/db'

export const config = {
  api: {
    bodyParser: false,
  },
}

export default async function Forward(req: NextApiRequest, res: NextApiResponse) {
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  const integrations = await db.objects.user_integrations.findUnique({
    where: {
      id: session.user.id
    }
  })
  if (integrations === null || !integrations.cavatica_api_key) throw new UnauthorizedError()
  const room_id = req.query.id
  const path = `/${((req.query.path ?? []) as string[]).join('/')}`
  const io = global.io
  if (io === undefined) throw new Error('Not Implemented')
  const request_id = randomUUID()
  const body = ['HEAD', 'GET'].includes(req.method ?? '') ? undefined : (await getRawBody(req)).toString('base64')
  const [socket, ..._] = await io.to(room_id).fetchSockets()
  const proxyReq = new Promise<{ status: number, body?: string, headers: Record<string, string> }>((resolve, reject) => {
    emitter.once(`http:recv:${request_id}`, (r) => resolve(r))
    setTimeout(() => reject(new TimeoutError()), 5000)
  })
  socket.emit('http:send', {
    id: request_id,
    path,
    headers: {
      'authorization': `Token ${integrations.cavatica_api_key}`,
      ...dict.fromIncomingHeaders(req.headers)
    },
    method: req.method,
    body: body,
  })
  try {
    const proxyRes = await proxyReq
    dict.items(proxyRes.headers).forEach(({ key, value }) => {
      res.setHeader(key, value)
    })
    res.status(proxyRes.status).send(proxyRes.body)
  } catch (e) {
    res.status(500).send(JSON.stringify(e))
  }
}
