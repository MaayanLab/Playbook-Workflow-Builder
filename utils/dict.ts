/**
 * Helpers to operate objects like lists (as can be done with python dicts)
 * We use {key, value} "Item" objects instead of tuples for improved compatibility with typescript.
 */

import type { IncomingHttpHeaders } from "http"

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
 * "Ordered" values
 */
export function keys<T extends {}>(dict: T): Array<keyof T> {
  return [...Object.keys(dict)] as Array<keyof T>
}

/**
 * "Ordered" values
 */
export function values<T extends {}>(dict: T): Array<T[keyof T]> {
  return keys(dict).map(key => dict[key])
}

/**
 * Python style object => list of tuples
 */
export function items<T extends {}>(dict: T) {
  return keys(dict).map((key) => ({ key, value: dict[key] }))
}

/**
 * Reduce dictionary by some filter, defaulting to ensuring the value evaluates to true
 */
export function filter<T extends {}, P extends (item: { key: keyof T, value: T[keyof T] }) => boolean>(dict: T, pred: P = (({ value }) => !!value) as P) {
  return init(items(dict).filter(pred))
}

/**
 * Check if a dictionary is empty
 */
export function isEmpty<T extends {}>(dict: T): boolean {
  return Object.keys(dict).length === 0
}

/**
 * "Sorted" keys
 */
export function sortedKeys<T extends {}>(dict: T): Array<keyof T> {
  const K = keys(dict)
  K.sort((a, b) => {
    const aStr = a.toString()
    const bStr = b.toString()
    if(aStr === bStr) {
      return 0
    } else if (aStr > bStr) {
      return 1
    } else {
      return -1
    }
  })
  return K
}

/**
 * "Sorted" values
 */
export function sortedValues<T extends {}>(dict: T): Array<T[keyof T]> {
  return sortedKeys(dict).map(key => dict[key])
}

/**
 * Python style object => list of tuples
 */
export function sortedItems<T extends {}>(dict: T) {
  return sortedKeys(dict).map((key) => ({ key, value: dict[key] }))
}

export function fromHeaders(headers: Headers) {
  const dict: Record<string, string> = {}
  headers.forEach((value, key) => {
    if ([
      'accept',
      'content-type',
      'authorization',
    ].includes(key.toLowerCase())) {
      dict[key] = value
    }
  })
  return dict
}

export function fromIncomingHeaders(headers: IncomingHttpHeaders) {
  return init(items(headers).filter(({ key }) => ([
    'accept' as keyof IncomingHttpHeaders,
    'content-type' as keyof IncomingHttpHeaders,
    'authorization' as keyof IncomingHttpHeaders,
  ].includes(key.toString().toLowerCase()))))
}
