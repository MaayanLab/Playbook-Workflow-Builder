export type MaybeArray<T> = T | T[]

export function ensureArray<T extends {}>(A?: T | T[]): T[] {
  if (typeof A === 'undefined') return []
  else if (typeof A === 'object' && Array.isArray(A)) {
    return A
  } else {
    return [A]
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
