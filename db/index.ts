/**
 * Database construction.
 *
 * DATABASE_URL selects the backend:
 *  - postgres://... or postgresql://...  postgres (production), migrate with `npm run db:migrate`
 *  - sqlite:<path>                      a sqlite database file, migrated automatically
 *  - unset                              an in-memory sqlite database, migrated automatically
 */
import { Kysely, PostgresDialect, SqliteDialect } from 'kysely'
import { Migrator } from 'kysely/migration'
import BetterSqlite3 from 'better-sqlite3'
import * as pg from 'pg'
import { randomUUID } from 'crypto'
import type { Database } from './schema'
import { migrationProvider } from './migrations'
import { wrapSqlite } from './sqlite'
import { type EventBus, LocalEventBus, PgEventBus } from './events'

export * from './schema'
export * from './constants'
export type { EventBus, Job } from './events'

export type DB = Kysely<Database>

export type Backend = {
  dialect: 'postgres' | 'sqlite'
  db: DB
  bus: EventBus
  /** the underlying pool when using postgres (used by the socket.io adapter) */
  pool?: pg.Pool
}

/** Generate a random id for a new record */
export const uuid = (): string => randomUUID()

export async function migrateToLatest(db: Kysely<any>) {
  const migrator = new Migrator({ db, provider: migrationProvider })
  const { error, results } = await migrator.migrateToLatest()
  for (const result of results ?? []) {
    if (result.status === 'Success') console.log(`applied migration ${result.migrationName}`)
    else if (result.status === 'Error') console.error(`failed to apply migration ${result.migrationName}`)
  }
  if (error) throw error
}

/**
 * Keep postgres values consistent with sqlite: timestamps (without time zone,
 *  stored as UTC) become ISO-8601 strings, and int8 becomes a javascript number.
 */
const pgTypes = {
  getTypeParser: ((oid: number, format?: any) => {
    if (oid === pg.types.builtins.TIMESTAMP) return (value: string) => new Date(value + 'Z').toISOString()
    if (oid === pg.types.builtins.INT8) return (value: string) => Number(value)
    return pg.types.getTypeParser(oid, format)
  }) as pg.CustomTypesConfig['getTypeParser'],
}

function createPostgres(connectionString: string): Backend {
  const pool = new pg.Pool({ connectionString, types: pgTypes })
  return {
    dialect: 'postgres',
    db: new Kysely<Database>({ dialect: new PostgresDialect({ pool }) }),
    bus: new PgEventBus(connectionString),
    pool,
  }
}

function createSqlite(filename: string): Backend {
  return {
    dialect: 'sqlite',
    db: new Kysely<Database>({
      dialect: new SqliteDialect({
        // the database is opened (and migrated) lazily on first use
        database: async () => {
          const raw = new BetterSqlite3(filename)
          if (filename !== ':memory:') raw.pragma('journal_mode = WAL')
          raw.pragma('foreign_keys = ON')
          const database = wrapSqlite(raw)
          await migrateToLatest(new Kysely<any>({ dialect: new SqliteDialect({ database }) }))
          return database
        },
      }),
    }),
    bus: new LocalEventBus(),
  }
}

export function createDatabase(url = ''): Backend {
  if (/^postgres(ql)?:\/\//.test(url)) return createPostgres(url)
  if (url.startsWith('sqlite:')) return createSqlite(url.slice('sqlite:'.length).replace(/^\/\//, '') || ':memory:')
  if (!url) return createSqlite(':memory:')
  throw new Error(`Unsupported DATABASE_URL, expected postgres://... or sqlite:<path>`)
}
