import React from 'react'

/**
 * Use localstorage as you would normal react state
 * 
 * const [myState, setMyState] = useLocalStorage('localStorageKey')
 * 
 * Only difference, it must be string|null
 *  & setMyState only takes a value not a callback
 */
export default function useLocalStorage(key: string) {
  const [value, _setValue] = React.useState(() => window.localStorage.getItem(key))
  const storageListener = React.useCallback((evt: StorageEvent | { key: string, newValue: string | null }) => {
    if (evt.key === key) {
      _setValue(() => evt.newValue)
    }
  }, [key])
  const setValue = React.useCallback((newValue: string | null) => {
    if (newValue !== null) {
      window.localStorage.setItem(key, newValue)
    } else {
      window.localStorage.removeItem(key)
    }
    storageListener({ key, newValue })
  }, [key])
  React.useEffect(() => {
    window.addEventListener('storage', storageListener)
    return () => {window.removeEventListener('storage', storageListener)}
  }, [key])
  return [value, setValue] as const
}