/**
 * A thin wrapper around better-sqlite3 for kysely's SqliteDialect which
 *  smooths over the differences between sqlite and postgres values:
 *  - booleans are bound as 0/1 and columns declared `boolean` are read back as booleans
 *  - Date parameters are bound as ISO-8601 strings
 */
import type BetterSqlite3 from 'better-sqlite3'
import type { SqliteDatabase, SqliteStatement } from 'kysely'

function bindable(value: unknown) {
  if (typeof value === 'boolean') return value ? 1 : 0
  if (value instanceof Date) return value.toISOString()
  return value
}

export function wrapSqlite(raw: BetterSqlite3.Database): SqliteDatabase {
  return {
    close: () => { raw.close() },
    prepare: (sql: string): SqliteStatement => {
      const stmt = raw.prepare(sql)
      const booleanColumns = stmt.reader
        ? stmt.columns().filter(col => col.type?.toLowerCase() === 'boolean').map(col => col.name)
        : []
      const decode = booleanColumns.length === 0
        ? (row: unknown) => row
        : (row: unknown) => {
          const record = row as Record<string, unknown>
          for (const col of booleanColumns) {
            if (record[col] !== null && record[col] !== undefined) record[col] = !!record[col]
          }
          return record
        }
      const bind = (parameters: ReadonlyArray<unknown>) => parameters.map(bindable)
      return {
        reader: stmt.reader,
        all: (parameters) => stmt.all(bind(parameters)).map(decode),
        run: (parameters) => stmt.run(bind(parameters)),
        iterate: function* (parameters) {
          for (const row of stmt.iterate(bind(parameters))) yield decode(row)
        },
      }
    },
  }
}
