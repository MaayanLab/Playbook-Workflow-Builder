import type { Socket } from 'socket.io'
import { EventEmitter } from 'stream'
import cache from '@/utils/global_cache'

export const emitter = cache('emitter', () => new EventEmitter())

export function onSocketCavatica(client: Socket) {
  client.on('cavatica:join', (id) => {
    console.debug(`${client.id} joined cavatica:${id}`)
    client.join(`cavatica:${id}`)
    emitter.emit(`join:${id}`)
    emitter.on(`close:${id}`, () => {
      console.debug(`${client.id} disconnected from ${id}`)
      client.leave(`cavatica:${id}`)
      client.send('cavatica:close')
    })
  })
  client.on('http:recv', ({ id, ...rest}) => {
    emitter.emit(`http:recv:${id}`, rest)
  })
}
