import type { Socket } from 'socket.io'
import fpprg from '@/app/fpprg'
import { TimeoutError } from '@/spec/error'
import { Resolved } from '@/core/FPPRG'

export function onSocketFPPRG(client: Socket) {
  client.on('fpprg:fpl', async (id: string) => {
    console.log(`requested fpl ${id}`)
    const fpl = await fpprg.getFPL(id)
    if (!fpl) {
      client.emit('fpprg:fpl', { id, error: 'fpl not found' })
    } else {
      console.log(`sent fpl ${id}`)
      client.emit('fpprg:fpl', { id, fpl: fpl.resolve().map(item => item.toJSON()) })
    }
  })
  client.on('fpprg:resolved', async (id: string) => {
    console.log(`requested resolved ${id}`)
    let resolved: Resolved | undefined
    while (resolved === undefined) {
      try {
        resolved = await fpprg.awaitResolved(id)
      } catch (e) {
        if (!(e instanceof TimeoutError)) {
          client.emit('fpprg:resolved', { id, error: 'An unexpected error occurred' } )
          console.error(e)
          throw e
        }
      }
      if (client.disconnected) return
    }
    console.log(`sent resolved ${id}`)
    client.emit('fpprg:resolved', resolved.toJSON())
  })
}
