import { Server } from 'socket.io'
import onSocket from '@/app/extensions/socket'

export const config = {
  api: {
    bodyParser: false
  }
}

export default function SocketHandler(req: any, res: any) {
  if (res.socket.server.io) {
    console.log('Socket is already running')
  } else {
    console.log('Socket is initializing')
    const io = new Server(res.socket.server)
    res.socket.server.io = io
    io.on('connection', onSocket)
  }
  res.end()
}
