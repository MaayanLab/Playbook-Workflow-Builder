export default function cache<T>(fn: (key: string) => T) {
  const _cache = new Map<string, T>()
  return (key: string, force = false) => {
    if (!_cache.has(key) || force) {
      _cache.set(key, fn(key))
    }
    return _cache.get(key) as T
  }
}