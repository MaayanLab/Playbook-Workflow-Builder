import http from 'http';
import { Server } from 'socket.io';
import db from '@/app/db'
import { createAdapter } from '@socket.io/postgres-adapter'
import { PgDatabase } from '@/utils/orm/pg'
import onSocket from '@/app/extensions/socket'

const server = http.createServer((req, res) => {
  // Handle HTTP requests if needed
});

const io = new Server(server, { transports: ['websocket'] });
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

server.listen(3005, '0.0.0.0', () => {
  console.log('WebSocket server listening on port 3005');
});
