export default function promise_cache<T>(fn: (key: string) => Promise<T>) {
  const _cache = new Map<string, Promise<T>>()
  return (key: string) => {
    if (!_cache.has(key)) {
      _cache.set(key, fn(key))
    }
    return _cache.get(key) as Promise<T>
  }
}