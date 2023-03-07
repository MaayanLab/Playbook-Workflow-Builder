import React from 'react'
import { BareFetcher, useSWRConfig, State } from 'swr'
import useAsyncEffect from 'use-async-effect'
import * as dict from '@/utils/dict'

/**
 * Apply useSWR to a list of keys, resolve all of them in parallel and return data as a Record<key, result>
 * We piggyback off of the SWR cache so that singular useSWR keys can also be cached.
 */
export default function useSWRMap<T>(keys: string[] | undefined, fetcher: BareFetcher<T>) {
  const {cache, mutate} = useSWRConfig()
  const [isLoading, setIsLoading] = React.useState<boolean>()
  const [data, setData] = React.useState<Record<string, T>>()
  const [error, setError] = React.useState<unknown>()
  useAsyncEffect(async (isMounted) => {
    if (!keys) return
    setIsLoading(() => true)
    try {
      const data = dict.init(
        await Promise.all(
          keys.map(async (key) => {
            let cachedValue = cache.get(key)
            if (cachedValue === undefined || !cachedValue.data) {
              await mutate(key, fetcher(key))
              cachedValue = cache.get(key) as State<any, any>
            }
            return { key, value: cachedValue.data }
          })
        )
      )
      if (!isMounted()) return
      setData(() => data)
      setError(() => undefined)
      setIsLoading(() => false)
    } catch (e) {
      if (!isMounted()) return
      setError(() => e)
      setIsLoading(() => false)
    }
  }, [(keys||[]).join('|')])
  return { data, error, isLoading }
}
