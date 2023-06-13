import emitter from '@/app/emitter'
import { Server } from 'socket.io'

export default function SocketHandler(req: any, res: any) {
  if (res.socket.server.io) {
    console.log('Socket is already running')
  } else {
    console.log('Socket is initializing')
    const io = new Server(res.socket.server)
    res.socket.server.io = io
    io.on('connection', (client) => {
      console.log(`${client} connected`)
      client.on('join', (room) => {
        client.join(room)
      })
      client.on('http:recv', ({ id, ...rest}) => {
        emitter.emit(`http:recv:${id}`, rest)
      })
    })
  }
  res.end()
}
