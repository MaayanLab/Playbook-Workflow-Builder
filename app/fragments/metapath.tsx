import React from 'react'
import type { inferRouterOutputs } from '@trpc/server'
import type { AppRouter } from '@/app/trpc/router'
import type KRG from '@/core/KRG'
import * as dict from '@/utils/dict'
import { useTRPCClient } from '@/app/trpc/client'
import { ResolutionWatcher } from '@/app/fragments/resolution-watcher'
import { useQuery } from '@tanstack/react-query'
import krg from '@/app/krg'
import { suggestionsQueryKey, withSuggestions } from '@/app/extensions/suggestions'

/** An FPL as received by the client (see FPL['toJSON']) */
export type Metapath = inferRouterOutputs<AppRouter>['fpprg']['fpl'][number]

type MetapathContextValue = {
  session_id?: string
  watcher: ResolutionWatcher
  krg: KRG
}

const MetapathContext = React.createContext<MetapathContextValue | undefined>(undefined)

export function MetapathProvider(props: React.PropsWithChildren<{ session_id?: string }>) {
  const trpc = useTRPCClient()
  const context = React.useMemo(() => ({ session_id: props.session_id }), [props.session_id])
  // One subscription transports every process watched on this page.
  const watcher = React.useMemo(() => new ResolutionWatcher(trpc, context), [trpc, context])
  const { data: suggestions } = useQuery({
    queryKey: suggestionsQueryKey(props.session_id),
    queryFn: () => trpc.catalog.suggestions.query(undefined, { context }),
    staleTime: 30_000,
    retry: false,
  })
  const sessionKrg = React.useMemo(() => withSuggestions(krg, suggestions ?? []), [suggestions])
  const value = React.useMemo(
    () => ({ session_id: props.session_id, watcher, krg: sessionKrg }),
    [props.session_id, watcher, sessionKrg],
  )
  return (
    <MetapathContext.Provider value={value}>
      {props.children}
    </MetapathContext.Provider>
  )
}

/** The metapath context if there is one (e.g. views also render in embeds without one) */
export function useOptionalMetapath() {
  return React.useContext(MetapathContext)
}

export function useMetapath() {
  const context = React.useContext(MetapathContext)
  if (!context) throw new Error('useMetapath must be used inside MetapathProvider')
  return context
}

/** The component registry and current database's user suggestions. */
export function useKRG() {
  return useMetapath().krg
}

export function useFPL(id?: string) {
  const { session_id } = useMetapath()
  const trpc = useTRPCClient()
  return useQuery<Metapath[]>({
    queryKey: ['fpprg-fpl', session_id, id],
    queryFn: async () => !id || id === 'start'
      ? []
      : await trpc.fpprg.fpl.query(id, { context: { session_id } }),
    staleTime: Infinity,
    placeholderData: previous => previous,
  })
}

export function useResolved(id?: string) {
  const { watcher } = useMetapath()
  const subscribe = React.useCallback((notify: () => void) => id ? watcher.watch(id, notify) : () => {}, [watcher, id])
  const getSnapshot = React.useCallback(() => id ? watcher.get(id) : undefined, [watcher, id])
  const data = React.useSyncExternalStore(subscribe, getSnapshot, () => undefined)
  return { data, mutate: () => { if (id) watcher.refresh(id) } }
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
  }, [resolved, krg, head.process.id, head.process.type])
  return { data: { output, outputNode }, error, isLoading, status, mutate }
}

/**
 * Retrieve inputs to this process from outputs of its inputs
 *  (requests are de-duplicated by the per-process cache in MetapathProvider)
 */
export function useMetapathInputs({ krg, head }: { krg: KRG, head: Metapath }) {
  const { watcher } = useMetapath()
  const keys = React.useMemo(() => dict.sortedKeys(head.process.inputs), [head.process.id])
  const ids = React.useMemo(() => keys.map(key => head.process.inputs[key].id), [keys, head.process.id])
  const subscribe = React.useCallback((notify: () => void) => watcher.watchMany(ids, notify), [watcher, ids])
  const version = React.useSyncExternalStore(subscribe, watcher.getVersion, () => 0)
  return React.useMemo(() => {
    const inputs: Record<string, any> = {}
    let error: Error | undefined
    keys.forEach((key, i) => {
      const value = watcher.get(ids[i])
      if (value?.type === 'error') error ??= new Error(value.error)
      const resolved = value?.type === 'resolved' ? value.data.data?.value : undefined
      if (resolved === undefined) return
      const [name, ...rest] = key.toString().split(':')
      if (rest.length > 0) {
        if (!(name in inputs)) inputs[name] = []
        inputs[name].push(resolved)
      } else {
        inputs[name] = resolved
      }
    })
    return { data: inputs, error }
  }, [watcher, version, keys, ids])
}
