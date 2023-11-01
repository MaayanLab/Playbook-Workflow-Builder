import React from 'react'
import type { FPL, Resolved } from '@/core/FPPRG'
import { Map } from 'immutable'
import type KRG from '@/core/KRG'
import * as dict from '@/utils/dict'
import { io } from 'socket.io-client'
import promise_cache from '@/utils/promise_cache'
import useSWR from 'swr'
import type { StatusUpdate } from '@/spec/metanode'

export type Metapath = ReturnType<FPL['toJSON']>
type ResolvedJSON = ReturnType<Resolved['toJSON']>

const MetapathContext = React.createContext({
  fetchFPL: (id: string, force?: boolean) => Promise.reject<Metapath[]>('Uninitialized'),
  fetchResolved: (id: string, force?: boolean) => Promise.reject<ResolvedJSON>('Uninitialized'),
})

export function MetapathProvider(props: React.PropsWithChildren<{ session_id?: string }>) {
  const fetchSocket = React.useCallback(promise_cache(async (session_id) => {
    await fetch(`${session_id ? `/api/${session_id}` : ''}/api/socket`)
    return io()
  }), [])
  const fetchFPL = React.useCallback(promise_cache(async (id) => {
    if (id === 'start') return []
    const socket = await fetchSocket(props.session_id ?? '')
    const ret = new Promise<Metapath[]>((resolve, reject) => {
      socket.once(`fpprg:fpl:${id}`, (ret: { error: any } | Metapath[]) => {
        if ('error' in ret) {
          reject(ret.error)
        } else {
          resolve(ret)
        }
      })
    })
    socket.emit('fpprg:fpl', id)
    return await ret
  }), [props.session_id])
  const fetchResolved = React.useCallback(promise_cache(async (id) => {
    const socket = await fetchSocket(props.session_id ?? '')
    const ret = new Promise<ResolvedJSON>((resolve, reject) => {
      socket.once(`fpprg:resolved:${id}`, (ret: { error: any } | ResolvedJSON) => {
        if ('error' in ret) {
          reject(ret.error)
        } else {
          resolve(ret)
        }
      })
    })
    const listener = (update: StatusUpdate) => {
      console.log(JSON.stringify({ [id]: update }))
    }
    socket.on(`fpprg:resolved:${id}:status`, listener)
    socket.emit('fpprg:resolved', id)
    try {
      return await ret
    } finally {
      socket.off(`fpprg:resolved:${id}:status`, listener)
    }
  }), [props.session_id])
  return (
    <MetapathContext.Provider value={{ fetchFPL, fetchResolved }}>
      {props.children}
    </MetapathContext.Provider>
  )
}

export function useMetapath() {
  return React.useContext(MetapathContext)
}

export function useFPL(id: string) {
  const metapath = useMetapath()
  const { mutate, ...swr } = useSWR(id, metapath.fetchFPL, { shouldRetryOnError: false })
  return { ...swr, mutate: async () => {await metapath.fetchFPL(id, true)} }
}

export function useResolved(id: string) {
  const metapath = useMetapath()
  const { mutate, ...swr } = useSWR(id, metapath.fetchResolved, { shouldRetryOnError: false })
  return { ...swr, mutate: async () => {await metapath.fetchResolved(id, true)} }
}

/**
 * Retrieve output from the API, decode and return
 */
export function useMetapathOutput({ krg, head }: { krg: KRG, head: Metapath }) {
  const { data, error, isLoading, mutate } = useResolved(head.process.id)
  const outputNode = React.useMemo(() => data?.data ? krg.getDataNode(data.data.type) : krg.getProcessNode(head.process.type).output, [data])
  const output = React.useMemo(() => data?.data === null ? null : data?.data.value, [data])
  return { data: { output, outputNode }, error, isLoading, mutate }
}

/**
 * Retrieve inputs to this process from outputs of its inputs
 *  We rely on SWR to help us de-duplicate these requests
 */
export function useMetapathInputs({ krg, head }: { krg: KRG, head: Metapath }) {
  const metapath = useMetapath()
  const [inputs, setInputs] = React.useState(Map<string, Map<string, any>>())
  React.useEffect(() => {
    setInputs(inputs => inputs.clear().set(head.process.id, Map()))
    dict.items(head.process.inputs)
      .forEach(({ key, value }) => {
        metapath.fetchResolved(value.id)
          .then(resolved => {
            setInputs(inputs => inputs.setIn([head.process.id, key], resolved.data))
          })
          .catch(error =>
            setInputs(inputs => inputs.setIn([head.process.id, key], { type: 'Error', value: error }))
          )
      })
  }, [head])
  const { error, isLoading } = React.useMemo(() =>
    inputs.get(head.id, Map())
      .valueSeq()
      .reduce(({ error, isLoading }, value: { type: string, data: any }) => ({
        error: value.type === 'Error' ? value.data : error,
        isLoading: isLoading || (value === undefined),
      }), { error: undefined as any, isLoading: true }),
  [head, inputs])
  return { data: inputs.get(head.process.id)?.toJS(), error, isLoading }
}
