import React from 'react'
import useSWR from 'swr'
import useSWRMutation from 'swr/mutation'
import * as dict from '@/utils/dict'
import fetcher, { fetcherGET } from '@/utils/next-rest-fetcher'
import type { APIRouteInterface } from "@/spec/api"

type FParamO<F> = F extends (param0: infer _P0, param1: infer _P1, opts: infer O) => infer _R ? O : never

export function useAPIQuery<Q extends {}, O>(route: APIRouteInterface<Q, O>, query: Q, { base = '', ...opts }: { base?: string } & FParamO<typeof useSWR<O>> = {} as { base?: string } & FParamO<typeof useSWR<O>>) {
  const key = React.useMemo(() => {
    let path = route.path
    const searchParams = new URLSearchParams()
    dict.items(query).forEach(({ key, value }) => {
      const m = (new RegExp(`\\[\\.*${key as string}\\]`)).exec(path)
      if (m !== null) path = path.replace(m[0], typeof value === 'string' ? value : JSON.stringify(value))
      else searchParams.append(key as string, typeof value === 'string' ? value : JSON.stringify(value))
    })
    const params = searchParams.toString()
    if (params) return path + '?' + params
    else return path
  }, [route.path, query])
  return useSWR(`${base ?? ''}${key}`, fetcherGET<O>, opts)
}

export function useAPIMutation<Q extends {}, O, B>(route: APIRouteInterface<Q, O, B>, query: Partial<Q> = {}, { base = '', ...opts }: { base?: string } & FParamO<typeof useSWRMutation<O, any, string, { query?: Partial<Q>, body?: B }>> = {} as { base?: string } & FParamO<typeof useSWRMutation<O, any, string, { query?: Partial<Q>, body?: B }>>) {
  return useSWRMutation(route.path, (_key: string, { arg = {} }: { arg?: { query?: Partial<Q>, body?: B } }) => {
    let path = route.path
    const searchParams = new URLSearchParams()
    dict.items({ ...query, ...(arg.query??{}) }).forEach(({ key, value }) => {
      const m = (new RegExp(`\\[\\.*${key as string}\\]`, 'g')).exec(path)
      if (m !== null) path = path.replace(m[0], typeof value === 'string' ? value : JSON.stringify(value))
      else searchParams.append(key as string, typeof value === 'string' ? value : JSON.stringify(value))
    })
    const params = searchParams.toString()
    if (params) path = path + '?' + params
    return fetcher<O>(`${base ?? ''}${path}`, { method: 'POST', body: typeof arg.body === 'undefined' ? undefined : arg.body instanceof FormData ? arg.body : JSON.stringify(arg.body) })
  }, opts)
}
