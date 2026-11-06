import React from 'react'
import type { inferRouterOutputs } from '@trpc/server'
import type { AppRouter } from '@/app/trpc/router'
import type KRG from '@/core/KRG'
import * as dict from '@/utils/dict'
import { useTRPCClient } from '@/app/trpc/client'
import cache from '@/utils/cache'
import { Readable, derived, readable } from '@/utils/store'
import { useReadable, usePromise } from '@/utils/hooks'
import type { ResolvedLifecycle } from '@/app/trpc/routers/fpprg'

/** An FPL as received by the client (see FPL['toJSON']) */
export type Metapath = inferRouterOutputs<AppRouter>['fpprg']['fpl'][number]

const MetapathContext = React.createContext({
  fetchFPL: (id: string) => Promise.reject<Metapath[]>('Uninitialized'),
  fetchResolved: (id: string) => Promise.reject<Readable<ResolvedLifecycle> & { refresh: () => void }>('Uninitialized'),
})

export function MetapathProvider(props: React.PropsWithChildren<{ session_id?: string }>) {
  const trpc = useTRPCClient()
  const context = React.useMemo(() => ({ session_id: props.session_id }), [props.session_id])
  const fetchFPL = React.useCallback(cache(async (id: string) => {
    if (id === 'start') return []
    return await trpc.fpprg.fpl.query(id, { context })
  }), [trpc, context])
  const fetchResolved = React.useCallback(cache(async (id: string) => {
    let refresh = () => {}
    const { subscribe } = readable<ResolvedLifecycle>(undefined, (set) => {
      let subscription: { unsubscribe: () => void } | undefined
      refresh = () => {
        subscription?.unsubscribe()
        subscription = trpc.fpprg.resolved.subscribe(id, {
          context,
          onData: (value) => { set(value) },
          onError: (error) => { set({ type: 'error', error: error.message }) },
        })
      }
      refresh()
      return () => { subscription?.unsubscribe() }
    })
    return { subscribe, refresh: () => refresh() }
  }), [trpc, context])
  return (
    <MetapathContext.Provider value={{ fetchFPL, fetchResolved }}>
      {props.children}
    </MetapathContext.Provider>
  )
}

export function useMetapath() {
  return React.useContext(MetapathContext)
}

export function useFPL(id?: string) {
  const metapath = useMetapath()
  const fpl = React.useMemo(() => id ? metapath.fetchFPL(id) : Promise.resolve([]), [id])
  return usePromise(fpl)
}

export function useResolved(id?: string) {
  const metapath = useMetapath()
  const resolved = React.useMemo(() => id ? metapath.fetchResolved(id) : Promise.resolve(undefined), [id])
  const { data: store } = usePromise(resolved)
  const data = useReadable(store)
  return { data, mutate: () => { store?.refresh() } }
}

/**
 * Retrieve output from the API, decode and return
 */
export function useMetapathOutput({ krg, head }: { krg: KRG, head: Metapath }) {
  const { data: resolved, mutate } = useResolved(head.process.id)
  const { output, outputNode, isLoading, status, error } = React.useMemo(() => {
    if (!resolved) {
      return {
        output: undefined,
        outputNode: krg.getProcessNode(head.process.type)?.output,
        status: undefined,
        isLoading: true,
        error: undefined,
      }
    }
    if (resolved.type !== 'resolved') {
      return {
        output: undefined,
        outputNode: krg.getProcessNode(head.process.type)?.output,
        isLoading: resolved.type === 'resolving',
        status: resolved.type === 'resolving' ? resolved.status : undefined,
        error: resolved.type === 'error' ? resolved.error : undefined,
      }
    }
    if (resolved.data.data !== null) {
      return {
        output: resolved.data.data.value,
        outputNode: krg.getDataNode(resolved.data.data.type),
        isLoading: false,
        status: undefined,
        error: undefined,
      }
    } else {
      return {
        output: null,
        outputNode: krg.getProcessNode(head.process.type).output,
        isLoading: false,
        status: undefined,
        error: undefined,
      }
    }
  }, [resolved])
  return { data: { output, outputNode }, error, isLoading, status, mutate }
}

/**
 * Retrieve inputs to this process from outputs of its inputs
 *  We rely on SWR to help us de-duplicate these requests
 */
export function useMetapathInputs({ krg, head }: { krg: KRG, head: Metapath }) {
  const metapath = useMetapath()
  const promise = React.useMemo(() => {
    const keys = dict.sortedKeys(head.process.inputs)
    return Promise.all(keys.map(key => metapath.fetchResolved(head.process.inputs[key].id)))
      .then(values => dict.init(values.map((value, i) => ({ key: keys[i], value }))))
  }, [metapath, head])
  const { data: stores, error } = usePromise(promise)
  const store = React.useMemo(() => {
    if (!stores) return
    const keys = dict.sortedKeys(stores)
    return derived(
      keys.map(key => stores[key]),
      (values) => {
        const inputs = {} as Record<string, any>
        values.forEach((value, i) => {
          const [k, ...rest] = keys[i].toString().split(':')
          const v = value.type === 'resolved' ? value.data.data?.value : undefined
          if (v === undefined) return
          if (rest.length > 0) {
            if (!(k in inputs)) inputs[k] = []
            inputs[k].push(v)
          } else {
            inputs[k] = v
          }
        })
        return inputs
      },
    )
  }, [stores, head.process.id])
  const data = useReadable(store)
  return { data, error }
}
