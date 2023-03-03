import useSWR from 'swr'
import useSWRMutation from 'swr/mutation'

type FParam0<F> = F extends (param0: infer P0, param1: infer P1, ...rest: infer A) => infer R ? P0 : never
type FParam1<F> = F extends (param0: infer P0, param1: infer P1, ...rest: infer A) => infer R ? P1 : never
type FParamR<F> = F extends (param0: infer P0, param1: infer P1, ...rest: infer A) => infer R ? A : never

export function useSWRPartial<T>(key: FParam0<typeof useSWR<T>>, fetcher: FParam1<typeof useSWR<T>>) {
  return (...rest: FParamR<typeof useSWR<T>>) => useSWR(key, fetcher, ...rest)
}

export function useSWRMutationPartial<T>(key: FParam0<typeof useSWRMutation<T>>, fetcher: FParam1<typeof useSWRMutation<T>>) {
  return (...rest: FParamR<typeof useSWRMutation<T>>) => useSWRMutation(key, fetcher, ...rest)
}
