import React, { PropsWithChildren } from 'react'
import type { FPL, Resolved } from '@/core/FPPRG'
import useSWRImmutable from 'swr/immutable'
import fetcher from '@/utils/next-rest-fetcher'
import useSWRMap from '@/utils/swr-map'
import { Error as ErrorComponent } from '@/components/core/error'
import type KRG from '@/core/KRG'
import * as dict from '@/utils/dict'
import { io, Socket } from 'socket.io-client'

export type Metapath = ReturnType<FPL['toJSON']>
type ResolvedJSON = ReturnType<Resolved['toJSON']>

const MetapathContext = React.createContext({
  fpl: undefined as Metapath[] | undefined,
  fplLookup: {} as Record<string, Metapath>,
  resolved: {} as Record<string, ResolvedJSON>,
  mutate: (processId: string) => {},
})

export function MetapathProvider({ children, session_id, id }: PropsWithChildren<{ session_id?: string, id: string }>) {
  const [fpls, setFPLs] = React.useState({} as Record<string, Metapath[]>)
  const [resolved, setResolved] = React.useState({} as Record<string, ResolvedJSON>)
  const socket = React.useRef<{
    socket: Socket,
    connected: boolean,
    queue: [string, any[]][],
  }>()
  const emit = React.useCallback((evt: string, ...args: any[]) => {
    if (!socket.current) return
    if (socket.current.connected) {
      socket.current.socket.emit(evt, ...args)
    } else {
      socket.current.queue.push([evt, args])
    }
  }, [socket])
  React.useEffect(() => {
    socket.current = { socket: io(), connected: false, queue: [] }
    fetch(`${session_id ? `/api/socket/${session_id}` : ''}/api/socket`).then(() => {
      if (!socket.current) return
      socket.current.socket.on('connect', () => {
        if (socket.current) {
          socket.current.connected = true
          while (socket.current.queue.length > 0 && socket.current.connected) {
            const [evt, args] = socket.current.queue.shift() as [string, any[]]
            socket.current.socket.emit(evt, ...args)
          }
        }
      })
      socket.current.socket.on('fpprg:fpl', (evt) => {
        if ('error' in evt) {
          console.error(evt)
        } else {
          setFPLs(() => ({ [evt.id]: evt.fpl }))
        }
      })
      socket.current.socket.on('fpprg:resolved', (evt) => {
        if ('error' in evt) {
          console.error(evt)
        } else {
          setResolved(current => ({ ...current, [evt.id]: evt }))
        }
      })
      socket.current.socket.on('disconnect', () => {
        if (socket.current) socket.current.connected = false
      })
    })
    return () => {
      if (socket.current) socket.current.socket.disconnect()
    }
  }, [socket, session_id])
  React.useEffect(() => {
    if (id === 'start') {
      setFPLs(() => ({ start: [] }))
    } else {
      emit('fpprg:fpl', id)
    }
  }, [id])
  const fpl = React.useMemo(() => fpls[id], [fpls, id])
  React.useEffect(() => {
    if (!fpl) return
    for (const el of fpl) {
      if (!(el.process.id in resolved)) {
        emit('fpprg:resolved', el.process.id)
      }
    }
  }, [fpl])
  const mutate = React.useCallback((processId: string) => {
    emit('fpprg:resolved', processId)
  }, [emit])
  const fplLookup = React.useMemo(() => fpl ? dict.init(fpl.map(el => ({ key: el.id, value: el }))) : {}, [fpl])
  return (
    <MetapathContext.Provider value={{ fpl, fplLookup, resolved, mutate }}>
      {children}
    </MetapathContext.Provider>
  )
}

export function useMetapath() {
  return React.useContext(MetapathContext)
}


/**
 * Retrieve output from the API, decode and return
 */
export function useMetapathOutput({ krg, head }: { krg: KRG, head: Metapath }) {
  const metapath = useMetapath()
  const rawOutput = React.useMemo(() => head.process.id in metapath.resolved ? metapath.resolved[head.process.id] : undefined, [metapath.resolved, head])
  const outputNode = React.useMemo(() => rawOutput?.data ? krg.getDataNode(rawOutput.data.type) : krg.getProcessNode(head.process.type).output, [rawOutput])
  const output = React.useMemo(() => rawOutput?.data ? rawOutput.data.value : undefined, [rawOutput])
  return { data: { output, outputNode }, error: undefined, isLoading: rawOutput?.data === undefined, mutate: () => metapath.mutate(head.process.id) }
}

/**
 * Retrieve inputs to this process from outputs of its inputs
 *  We rely on SWR to help us de-duplicate these requests
 */
export function useMetapathInputs({ krg, head }: { krg: KRG, head: Metapath }) {
  const metapath = useMetapath()
  const fpl = React.useMemo(() => metapath.fplLookup[head.id], [metapath.fplLookup, head])
  const inputs = React.useMemo(() => fpl ? dict.init(dict.items(fpl.process.inputs).map(({ key, value }) => ({ key, value: metapath.resolved[value.id] ? metapath.resolved[value.id].data?.value : undefined }))) : undefined, [metapath.resolved, fpl])
  const isLoading = React.useMemo(() => inputs ? dict.values(inputs).some(value => value === undefined) : true, [inputs])
  return { data: inputs, error: undefined, isLoading }
}
