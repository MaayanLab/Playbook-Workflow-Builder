import db from '@/app/db'
import { Server } from 'socket.io'
import { createAdapter } from '@socket.io/postgres-adapter'
import { PgDatabase } from '@/utils/orm/pg'
import onSocket from '@/app/extensions/socket'

export default function SocketHandler(req: any, res: any) {
  if (res.socket.server.io) {
    console.log('Socket is already running')
  } else {
    console.log('Socket is initializing')
    const io = new Server(res.socket.server, { transports: ['websocket'] })
    if ('pool' in db) {
      // if there would be multiple instances of the ui, this will ensure
      //  messages can be shared across processes
      (db as PgDatabase).pool.query(`
        CREATE TABLE IF NOT EXISTS socket_io_attachments (
          id          bigserial UNIQUE,
          created_at  timestamptz DEFAULT NOW(),
          payload     bytea
        );
      `)
      io.adapter(createAdapter((db as PgDatabase).pool))
    }
    res.socket.server.io = io
    io.on('connection', onSocket)
  }
  res.end()
}
