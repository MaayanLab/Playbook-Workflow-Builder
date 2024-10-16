import type { MaybeArray } from '@/utils/types'

export function ensureArray<T extends {}>(A?: MaybeArray<T>): Array<T> {
  if (typeof A === 'undefined') return []
  else if (typeof A === 'object' && Array.isArray(A)) {
    return A
  } else {
    return [A]
  }
}

export function ensureOne<T>(A: T | T[]): T {
  if (typeof A === 'object' && Array.isArray(A)) {
    return A[0]
  } else {
    return A
  }
}

export function arange(n: number) {
  let R = []
  for (let i = 0; i < n; i++) {
    R.push(i)
  }
  return R
}

export function unique<T>(array: T[]) {
  const set = new Set<T>()
  array.forEach(element => set.add(element))
  return Array.from(set)
}

export function intersection<T>(As: T[], Bs: T[]) {
  const A = new Set<T>()
  const B = new Set<T>()
  As.forEach(a => A.add(a))
  Bs.forEach(b => B.add(b))
  return Array.from(A).filter(a => B.has(a))
}

export function any<T>(L: T[]) {
  return L.some(el => !!el)
}

export function all<T>(L: T[]) {
  return !L.some(el => !el)
}

export function chunked<T>(L: T[], chunk_size = 100) {
  const chunked: T[][] = []
  let current_chunk: T[] = []
  for (const el of L) {
    current_chunk.push(el)
    if (current_chunk.length === chunk_size) {
      chunked.push(current_chunk)
      current_chunk = []
    }
  }
  if (current_chunk.length > 0) chunked.push(current_chunk)
  return chunked
}
