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
      console.debug(`${client.id} connected`)
      client.on('join', (id) => {
        console.debug(`${client.id} joined ${id}`)
        client.join(id)
        emitter.emit(`join:${id}`)
        emitter.on(`close:${id}`, () => {
          console.debug(`${client.id} disconnected from ${id}`)
          client.leave(id)
          client.send('close')
        })
      })
      client.on('http:recv', ({ id, ...rest}) => {
        emitter.emit(`http:recv:${id}`, rest)
      })
    })
  }
  res.end()
}
