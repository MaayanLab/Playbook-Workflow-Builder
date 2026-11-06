/**
 * Migrations are registered explicitly (rather than read from disk)
 *  so that they work when bundled by next.js.
 *
 * Add new migrations here with a lexicographically increasing key.
 */
import type { Migration, MigrationProvider } from 'kysely/migration'
import * as m0000_baseline from './0000_baseline'

const migrations: Record<string, Migration> = {
  '0000_baseline': m0000_baseline,
}

export const migrationProvider: MigrationProvider = {
  getMigrations: async () => migrations,
}
