import type { Socket, Server } from 'socket.io'
import { EventEmitter } from 'stream'
import cache from '@/utils/global_cache'

export const emitter = cache('emitter', () => new EventEmitter())

export function onSocketCavatica(server: Server, client: Socket) {
  client.on('join', async (id, ack) => {
    // console.debug(`[${id}]: client ${client.id} joined`)
    await client.join(`client:${id}`)
    if (ack) ack()
  })
  client.on('leave', async (id, ack) => {
    // console.debug(`[${id}]: client ${client.id} left`)
    await client.leave(`client:${id}`)
    if (ack) ack()
  })
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
  client.on('http:recv', ({ id, ...rest}) => {
    emitter.emit(`http:recv:${id}`, rest)
  })
  client.onAny(async (evt, ...args) => {
    const evt_m = /^ws:(.+)$/.exec(evt)
    if (evt_m === null) return
    client.rooms.forEach(room => {
      const room_m = /^(client|worker):(.+)$/.exec(room)
      if (room_m === null) return
      if (room_m[1] === 'client') {
        // console.debug(`[${room_m[1]}]: forward ${evt} to workers`)
        server.to(`worker:${room_m[2]}`).emit(evt, ...args)
      } else if (room_m[1] === 'worker') {
        // console.debug(`[${room_m[1]}]: forward ${evt} to clients`)
        server.to(`client:${room_m[2]}`).emit(evt, ...args)
      }
    })
  })
}
