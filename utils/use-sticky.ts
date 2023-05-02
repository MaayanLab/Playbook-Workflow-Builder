import React from 'react'
import useSWRImmutable from 'swr/immutable'

export function useSticky(value: any) {
  const ref = React.useRef()
  if (ref !== undefined && value !== undefined) ref.current = value
  return ref.current
}

export function useSWRImmutableSticky<T = any>(key?: string | (() => string|undefined)) {
  const swr = useSWRImmutable<T>(key)
  const stickyData = useSticky(swr.data)
  return { ...swr, data: swr.data === undefined ? stickyData : swr.data, isLoading: swr.data === undefined }
}
