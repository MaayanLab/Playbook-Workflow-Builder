export type MaybeArray<T> = T | T[]
export type ValuesOfMaybeArray<T> = T extends (infer T_)[] ? T_ : T
export type ValuesOf<T> = T[keyof T]

/**
 * This helper is used extensively: it asserts that type parameter K exists as a key in T and returns it. This is
 *  a useful trick for storing lots of type information in T, organized by K
 */
export type ExtractKey<T, K> = T extends {[k in keyof T]: T[k]} ? K extends keyof T ? T[K] : never : never

/**
 * This helper lets us add inline type bounds
 */
export type Ensure<T, B> = T extends B ? T : never
