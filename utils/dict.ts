/**
 * Helpers to operate objects like lists (as can be done with python dicts)
 * We use {key, value} "Item" objects instead of tuples for improved compatibility with typescript.
 */

/**
 * Python style object => list of tuples
 */
export function items<T extends {}>(dict: T) {
  const K = keys(dict)
  return K.map((key) => ({ key, value: dict[key] }))
}

/**
 * Python style list of tuples => object
 */
export function init<K extends string | number | symbol, V>(items: { key: K, value: V}[]) {
  const output = {} as Record<K, V>
  for (const { key, value } of items) {
    output[key] = value
  }
  return output
}

/**
 * "Stable" keys
 */
export function keys<T extends {}>(dict: T): Array<keyof T> {
  const keys = Object.keys(dict) as Array<keyof T>
  keys.sort()
  return keys
}

/**
 * "Stable" values
 */
export function values<T extends {}>(dict: T): Array<T[keyof T]> {
  return keys(dict).map(key => dict[key])
}

/**
 * Reduce dictionary by some filter, defaulting to ensuring the value evaluates to true
 */
export function filter<T extends {}, P extends (item: { key: keyof T, value: T[keyof T] }) => boolean>(dict: T, pred: P = (({ value }) => !!value) as P) {
  return init(items(dict).filter(pred))
}