/**
 * Helpers to operate objects like lists (as can be done with python dicts)
 * We use {key, value} "Item" objects instead of tuples for improved compatibility with typescript.
 */
type Item<T = unknown> = { key: string|number|symbol, value: T }

/**
 * Python style object => list of tuples
 */
export function items<T>(dict: Record<string|number|symbol, T>) {
  const K = Object.keys(dict)
  K.sort() // for stability
  return K.map((key) => ({ key, value: dict[key] } as Item<T>))
}

/**
 * Python style list of tuples => object
 */
export function init<T>(items: Item<T>[]) {
  const output: Record<string|number|symbol, T> = {}
  for (const { key, value } of items) {
    output[key] = value
  }
  return output
}

/**
 * "Stable" keys
 */
export function keys<T>(dict: {[K in keyof T]: unknown}): Array<keyof T> {
  const keys = Object.keys(dict)
  keys.sort()
  return keys as Array<keyof T>
}

/**
 * "Stable" values
 */
export function values<T extends {}>(dict: T): Array<T[keyof T]> {
  return keys(dict).map(key => dict[key])
}
