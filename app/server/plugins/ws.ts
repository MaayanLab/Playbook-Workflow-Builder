import http from 'http';
import { Server } from 'socket.io';
import db from '@/app/db'
import { createAdapter } from '@socket.io/postgres-adapter'
import { PgDatabase } from '@/utils/orm/pg'
import onSocket from '@/app/extensions/socket'
import '@/utils/global_cache'
declare global {
  var io: Server | undefined
}

export default async function plugin(server: http.Server, opts: {}) {
  console.log('Starting Socket.io server...')
  const io = global.io = new Server(server, { transports: ['websocket'] });
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
  io.on('connection', (socket) => onSocket(io, socket));
}
