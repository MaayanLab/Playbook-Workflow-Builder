import type { Socket, Server } from 'socket.io'
import { EventEmitter } from 'stream'
import cache from '@/utils/global_cache'

export const emitter = cache('emitter', () => new EventEmitter())

export function onSocketCavatica(server: Server, client: Socket) {
  client.on('cavatica:join', async (id) => {
    console.debug(`${client.id} joined cavatica:${id}`)
    await client.join(`cavatica:${id}`)
    emitter.emit(`join:${id}`)
    emitter.on(`close:${id}`, async () => {
      console.debug(`${client.id} disconnected from ${id}`)
      await client.leave(`cavatica:${id}`)
      client.send('cavatica:close')
    })
  })
  client.on('http:recv', ({ id, ...rest}) => {
    emitter.emit(`http:recv:${id}`, rest)
  })
  client.onAny(async (evt, ...args) => {
    const m = /^ws:([^:]+):(.+)$/.exec(evt)
    if (m === null) return
    const sockets = await server.in(m[1]).fetchSockets()
    sockets.filter(socket => socket.id !== client.id).forEach(socket => {
      socket.emit(evt, ...args)
    })
  })
}
