/**
 * CAVATICA workers connect to the public server over socket.io and join `worker:${session_id}`.
 *  - http requests to /api/socket/[session_id]/... are forwarded with `http:send`/`http:recv`
 *  - websocket connections to the tRPC endpoint with ?session_id= are tunneled with `tunnel:*`
 */
import type { Socket, Server } from 'socket.io'
import type { IncomingMessage } from 'http'
import type { Duplex } from 'stream'
import { EventEmitter } from 'events'
import { randomUUID } from 'crypto'
import { WebSocketServer, type WebSocket } from 'ws'
import db from '@/app/db'
import cache from '@/utils/global_cache'
import { getSessionFromRequest } from '@/app/trpc/context'

export const emitter = cache('emitter', () => new EventEmitter())

const tunnelWss = cache('tunnel-wss', () => new WebSocketServer({ noServer: true }))
const tunnels = cache('tunnels', () => new Map<string, { session_id: string, ws: WebSocket }>())

function reject(socket: Duplex, status: string) {
  socket.write(`HTTP/1.1 ${status}\r\nConnection: close\r\n\r\n`)
  socket.destroy()
}

/**
 * Accept a browser websocket and relay its frames to the worker for `session_id`
 */
export async function tunnelUpgrade(server: Server, session_id: string, req: IncomingMessage, socket: Duplex, head: Buffer) {
  const session = await getSessionFromRequest(req)
  if (!session?.user) return reject(socket, '401 Unauthorized')
  const integrations = await db.selectFrom('user_integrations').select('cavatica_api_key').where('id', '=', session.user.id).executeTakeFirst()
  if (!integrations?.cavatica_api_key) return reject(socket, '401 Unauthorized')
  tunnelWss.handleUpgrade(req, socket, head, (ws) => {
    const conn = randomUUID()
    const worker = server.to(`worker:${session_id}`)
    tunnels.set(conn, { session_id, ws })
    worker.emit('tunnel:open', { conn })
    ws.on('message', (data) => { worker.emit('tunnel:send', { conn, data: data.toString() }) })
    ws.on('close', () => {
      if (tunnels.delete(conn)) worker.emit('tunnel:close', { conn })
    })
  })
}

export function onSocketCavatica(server: Server, client: Socket) {
  /** only accept tunnel messages from the worker serving that tunnel's session */
  const tunnelFor = (conn: unknown) => {
    const tunnel = typeof conn === 'string' ? tunnels.get(conn) : undefined
    if (tunnel && client.rooms.has(`worker:${tunnel.session_id}`)) return tunnel
  }
  client.on('worker:join', async (id, ack) => {
    // console.debug(`[${id}]: worker ${client.id} joined`)
    await client.join(`worker:${id}`)
    emitter.emit(`join:${id}`)
    emitter.on(`close:${id}`, async () => {
      // console.debug(`[${id}]: worker ${client.id} disconnected`)
      await client.leave(`worker:${id}`)
      client.emit('cavatica:close')
    })
    if (ack) ack()
  })
  client.on('http:recv', ({ id, ...rest }) => {
    emitter.emit(`http:recv:${id}`, rest)
  })
  client.on('tunnel:recv', ({ conn, data }) => {
    tunnelFor(conn)?.ws.send(data)
  })
  client.on('tunnel:close', ({ conn }) => {
    const tunnel = tunnelFor(conn)
    if (tunnel) {
      tunnels.delete(conn)
      tunnel.ws.close()
    }
  })
  client.on('disconnecting', () => {
    // close tunnels to this worker
    for (const [conn, tunnel] of tunnels) {
      if (client.rooms.has(`worker:${tunnel.session_id}`)) {
        tunnels.delete(conn)
        tunnel.ws.close()
      }
    }
  })
}
