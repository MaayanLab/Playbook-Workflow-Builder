import React from 'react'
import { Readable } from './store'

export function usePromise<T>(promise: Promise<T>) {
  const [value, setValue] = React.useState<{ data?: T, error?: any }>({})
  React.useEffect(() => {
    setValue(() => ({}))
    promise
      .then(data => setValue(() => ({ data })))
      .catch(error => setValue(() => ({ error })))
    return () => {setValue(() => ({}))}
  }, [promise])
  return value
}

export function useReadable<T>(store?: Readable<T>) {
  const [value, setValue] = React.useState<T>()
  React.useEffect(() => {
    if (!store) setValue(undefined)
    else return store.subscribe(value => {setValue(() => value)})
  }, [store])
  return value
}
