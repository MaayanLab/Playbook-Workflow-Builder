export type MaybeArray<T> = T | T[]
export type ValuesOfMaybeArray<T> = T extends (infer T_)[] ? T_ : T
export type ValuesOf<T> = T[keyof T]