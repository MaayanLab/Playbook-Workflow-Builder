import type { Socket, Server } from 'socket.io'
import fpprg from '@/app/fpprg'
import { TimeoutError } from '@/spec/error'
import { Resolved } from '@/core/FPPRG'
import db from '@/app/db'
import { StatusUpdate } from '@/spec/metanode'

export type ResolvedJSON = ReturnType<Resolved['toJSON']>
export type ResolvedLifecycle = {
  type: 'error',
  error: string,
} | {
  type: 'resolving',
  percent: number,
  status: string,
} | {
  type: 'resolved',
  data: ResolvedJSON
}

export function onSocketFPPRG(server: Server, client: Socket) {
  client.on('fpprg:fpl', async (id: string) => {
    const fpl = await fpprg.getFPL(id)
    if (!fpl) {
      client.emit(`fpprg:fpl:${id}`, { error: 'fpl not found' })
    } else {
      client.emit(`fpprg:fpl:${id}`, fpl.resolve().map(item => item.toJSON()))
    }
  })
  client.on('fpprg:resolved', async (id: string) => {
    const status = {
      type: 'resolving' as const,
      percent: 0,
      status: '',
    }
    // forward partial resolution update to the client
    const unsub = db.listen(`distributed:resolved:${id}:status`, (rawUpdate: string) => {
      const update: StatusUpdate = JSON.parse(rawUpdate)
      if (update.type === 'progress') status.percent = update.percent
      else if (update.type === 'info') status.status += update.message
      client.emit(`fpprg:resolved:${id}`, status)
    })
    // wait for it to be resolved
    let resolved: Resolved | undefined
    while (resolved === undefined) {
      try {
        resolved = await fpprg.awaitResolved(id)
      } catch (e) {
        if (!TimeoutError.isinstance(e)) {
          client.emit(`fpprg:resolved:${id}`, { type: 'error', error: (e as Error).toString() } )
          console.error(e)
          unsub()
          throw e
        }
      }
      if (client.disconnected) {
        unsub()
        return
      }
    }
    unsub()
    // send the resolved instance to the requesting user
    client.emit(`fpprg:resolved:${resolved.id}`, {
      type: 'resolved',
      data: resolved.toJSON(),
    })
  })
}
