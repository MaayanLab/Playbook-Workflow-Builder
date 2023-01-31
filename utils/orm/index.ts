import * as pg from 'pg'
import { PgTable } from "./pg"
import { MemoryTable } from "./memory"
import { DbTable } from './common'
import * as dict from '@/utils/dict'
import type { TypedSchema } from "@/spec/sql"
import { Codec } from '@/spec/codec'

export type DbOptions<T extends {}> = {
  connectionString?: string,
  schema: {[K in keyof T]: TypedSchema<T[K]>},
}

export type DbDatabase<T extends {[K in keyof T]: T[K]}> = {[K in keyof T]: DbTable<T[K] extends { id: Codec<string, string> } ? T[K] : never>}

export default function create_database<T extends {[K in keyof T]: T[K]}>(options: DbOptions<T>): DbDatabase<T> {
  if (options.connectionString) {
    const pool = new pg.Pool({ connectionString: options.connectionString })
    return dict.init(dict.items(options.schema).map(({ key, value }) => ({ key, value: new PgTable(value as any, pool) }))) as DbDatabase<T>
  } else {
    const memory: any = {}
    return dict.init(dict.items(options.schema).map(({ key, value }) => ({ key, value: new MemoryTable(value as any, memory) }))) as DbDatabase<T>
  }
}
