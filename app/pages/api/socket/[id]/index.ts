import type { Server } from "socket.io"
import type { NextApiRequest, NextApiResponse } from "next"
import { randomUUID } from "crypto"

export default function Forward(req: NextApiRequest, res: NextApiResponse) {
  const room_id = req.query.id as string
  const io = (res as any).socket.server.io as Server
  const request_id = randomUUID()
  const backend = io.in(room_id)
  io.once(`http:recv:${request_id}`, (r) => {
    res.status(r.status).send(r.body)
  })
  backend.emit('http:send', { id: request_id, path: req.method, body: req.body })
}
