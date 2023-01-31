declare global {
  var cache: any
}

export default function cache<T>(key: string, value: T): T {
  if (!global.cache) global.cache = {}
  if (!(key in global.cache)) global.cache[key] = value
  return global.cache[key] as T
}
