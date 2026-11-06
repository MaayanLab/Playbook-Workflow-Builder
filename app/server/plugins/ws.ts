import http from 'http'
import { Server } from 'socket.io'
import { WebSocketServer } from 'ws'
import { applyWSSHandler } from '@trpc/server/adapters/ws'
import { createAdapter } from '@socket.io/postgres-adapter'
import { pool } from '@/app/db'
import onSocket from '@/app/extensions/socket'
import { tunnelUpgrade } from '@/app/extensions/socket/cavatica'
import { appRouter } from '@/app/trpc/router'
import { createWSContext } from '@/app/trpc/context'
import { TRPC_WS_PATH } from '@/app/trpc/constants'
import { claimUpgradePath } from '../upgrade'
import '@/utils/global_cache'
declare global {
  var io: Server | undefined
}

export default async function plugin(server: http.Server, opts: {}) {
  // socket.io is used for the tunnel between this server & CAVATICA workers
  console.log('Starting Socket.io server...')
  const io = global.io = new Server(server, {
    transports: ['websocket'],
    // leave other upgrade requests (i.e. tRPC) alone
    destroyUpgrade: false,
  })
  if (pool) {
    // if there would be multiple instances of the ui, this will ensure
    //  messages can be shared across processes
    await pool.query(`
      CREATE TABLE IF NOT EXISTS socket_io_attachments (
        id          bigserial UNIQUE,
        created_at  timestamptz DEFAULT NOW(),
        payload     bytea
      );
    `)
    io.adapter(createAdapter(pool))
  }
  io.on('connection', (socket) => onSocket(io, socket))
  claimUpgradePath('/socket.io/')

  // tRPC subscriptions
  console.log('Starting tRPC websocket server...')
  const wss = new WebSocketServer({ noServer: true })
  const trpcHandler = applyWSSHandler({
    wss,
    router: appRouter,
    createContext: createWSContext,
    keepAlive: { enabled: true, pingMs: 30_000, pongWaitMs: 5_000 },
  })
  claimUpgradePath(TRPC_WS_PATH)
  server.on('upgrade', (req, socket, head) => {
    const url = new URL(req.url ?? '/', 'http://localhost')
    if (url.pathname !== TRPC_WS_PATH) return
    const session_id = url.searchParams.get('session_id')
    if (session_id) {
      // forward to the CAVATICA worker handling this session
      tunnelUpgrade(io, session_id, req, socket, head).catch((error) => {
        console.error(error)
        socket.destroy()
      })
    } else {
      wss.handleUpgrade(req, socket, head, (ws) => { wss.emit('connection', ws, req) })
    }
  })
  process.on('SIGTERM', () => {
    trpcHandler.broadcastReconnectNotification()
    wss.close()
  })
}
