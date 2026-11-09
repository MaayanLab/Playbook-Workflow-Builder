/**
 * Copy-on-write helpers for using a native Map as React state:
 *  they return the same map when nothing changes so React can bail out of re-rendering.
 */
export function mapSet<K, V>(map: ReadonlyMap<K, V>, key: K, value: V): ReadonlyMap<K, V> {
  if (map.has(key) && map.get(key) === value) return map
  return new Map(map).set(key, value)
}

export function mapDelete<K, V>(map: ReadonlyMap<K, V>, key: K): ReadonlyMap<K, V> {
  if (!map.has(key)) return map
  const copy = new Map(map)
  copy.delete(key)
  return copy
}

export function mapValues<K, V, U>(map: ReadonlyMap<K, V>, fn: (value: V, key: K) => U): Map<K, U> {
  return new Map([...map].map(([key, value]) => [key, fn(value, key)]))
}
