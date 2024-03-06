import { PgDatabase, PgTable } from "./pg"
import { MemoryDatabase, MemoryTable } from "./memory"
import * as dict from '@/utils/dict'
import type { TypedSchema } from "@/spec/sql"
import { DbDatabase, DbTables } from "./common"

export type DbOptions<T extends {}> = {
  connectionString?: string,
  schema: {[K in keyof T]: TypedSchema<T[K]>},
}

export type Db<T extends {}> = Omit<DbDatabase, 'objects'> & { objects: DbTables<T> }

export default function create_database<T extends {[K in keyof T]: T[K]}>(options: DbOptions<T>): Db<T> {
  if (options.connectionString) {
    const db = new PgDatabase(options.connectionString)
    const objects = dict.init(dict.items(options.schema).map(({ key, value }) => ({ key, value: new PgTable(value, db) })))
    db.objects = objects
    return { ...db, objects }
  } else {
    const db = new MemoryDatabase()
    const objects = dict.init(dict.items(options.schema).map(({ key, value }) => ({ key, value: new MemoryTable(value, db) })))
    db.objects = objects
    return { ...db, objects }
  }
}
