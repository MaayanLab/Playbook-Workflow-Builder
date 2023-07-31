import type { Decoded } from '@/spec/codec'
import type { TypedSchema, TypedSchemaRecord } from '@/spec/sql'

export type Data<T> = {[K in keyof T]: Decoded<T[K]>}
export type PartialData<T> = Partial<{[K in keyof T]?: Decoded<T[K]>}>
export type Where<T> = Partial<{[K in keyof T]?: Decoded<T[K]>}>
export type WhereMany<T> = Partial<{[K in keyof T]?: { in: Decoded<T[K]>[] } | Decoded<T[K]> }>
export type OrderBy<T> = Partial<{[K in keyof T]?: 'asc' | 'desc'}>
export type Find<T> = {
  select?: Partial<{[K in keyof T]?: boolean}>
  where: Where<T>
}
export type FindMany<T> = {
  select?: Partial<{[K in keyof T]?: boolean}>
  where?: WhereMany<T>
  orderBy?: OrderBy<T>
  skip?: number
  take?: number
}
export type Create<T> = {
  data: PartialData<T>
}
export type Update<T> = {
  where: Where<T>
  data: PartialData<T>
}
export type Upsert<T> = {
  where: Where<T>
  update?: PartialData<T>
  create: PartialData<T>
}
export type Delete<T> = {
  where: Where<T>
}

export interface DbTable<T extends {}> {
  create: (create: Create<T>) => Promise<TypedSchemaRecord<TypedSchema<T>>>
  findUnique: (find: Find<T>) => Promise<TypedSchemaRecord<TypedSchema<T>> | null>
  findMany: (find?: FindMany<T>) => Promise<Array<TypedSchemaRecord<TypedSchema<T>>>>
  update: (update: Update<T>) => Promise<TypedSchemaRecord<TypedSchema<T>> | null>
  upsert: (upsert: Upsert<T>) => Promise<TypedSchemaRecord<TypedSchema<T>>>
  delete: (delete_: Delete<T>) => Promise<TypedSchemaRecord<TypedSchema<T>> | null>
}

export interface DbDatabase {
  objects: any

  // subscriber
  listen: (cb: (evt: string, data: unknown) => void) => () => void

  // boss
  send: (queue: string, work: unknown) => Promise<void>
  work: (queue: string, opts: unknown, cb: (work: unknown) => Promise<void>) => Promise<() => void>
}

export type DbTables<T> = {[K in keyof T]: DbTable<T[K] extends {} ? T[K] : never>}

export type DbOptions<T extends {}> = {
  schema: {[K in keyof T]: TypedSchema<T[K]>},
}