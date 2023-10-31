import type { Socket } from 'socket.io'
import fpprg from '@/app/fpprg'
import { TimeoutError } from '@/spec/error'
import { Resolved } from '@/core/FPPRG'

export function onSocketFPPRG(client: Socket) {
  client.on('fpprg:fpl', async (id: string) => {
    const fpl = await fpprg.getFPL(id)
    if (!fpl) {
      client.emit(`fpprg:fpl:${id}`, { error: 'fpl not found' })
    } else {
      client.emit(`fpprg:fpl:${id}`, fpl.resolve().map(item => item.toJSON()))
    }
  })
  client.on('fpprg:resolved', async (id: string) => {
    let resolved: Resolved | undefined
    while (resolved === undefined) {
      try {
        resolved = await fpprg.awaitResolved(id)
      } catch (e) {
        if (!(e instanceof TimeoutError)) {
          client.emit(`fpprg:resolved:${id}`, { error: 'An unexpected error occurred' } )
          console.error(e)
          throw e
        }
      }
      if (client.disconnected) return
    }
    client.emit(`fpprg:resolved:${resolved.id}`, resolved.toJSON())
  })
}
