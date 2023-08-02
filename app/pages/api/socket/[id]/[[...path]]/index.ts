import type { Server } from "socket.io"
import type { NextApiRequest, NextApiResponse } from "next"
import { randomUUID } from "crypto"
import getRawBody from 'raw-body'
import { TimeoutError } from "@/spec/error"
import emitter from "@/app/emitter"

export const config = {
  api: {
    bodyParser: false,
  },
}

export default async function Forward(req: NextApiRequest, res: NextApiResponse) {
  const room_id = req.query.id as string
  const path = `/${((req.query.path ?? []) as string[]).join('/')}`
  const io = (res as any).socket.server.io as Server
  const request_id = randomUUID()
  const body = (await getRawBody(req)).toString()
  const [socket, ..._] = await io.to(room_id).fetchSockets()
  const proxyReq = new Promise<{ status: number, body?: string, headers: any }>((resolve, reject) => {
    emitter.once(`http:recv:${request_id}`, (r) => resolve(r))
    setTimeout(() => reject(new TimeoutError()), 5000)
  })
  socket.emit('http:send', { id: request_id, path, method: req.method, body: body ? body : undefined })
  try {
    const proxyRes = await proxyReq
    res.status(proxyRes.status).send(proxyRes.body)
  } catch (e) {
    res.status(500).send(JSON.stringify(e))
  }
}
