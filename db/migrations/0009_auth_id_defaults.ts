/**
 * With `advanced.database.generateId: 'uuid'` (app/auth.ts), Better Auth leaves ids to the database
 *  where it supports uuids (postgres), expecting a default, as its own schema has:
 *  `id ... DEFAULT gen_random_uuid()`. The auth tables created by 0006 & 0008 had no default, so
 *  inserting into them failed on postgres (sqlite doesn't support uuids: Better Auth generates them).
 */
import { type Kysely, PostgresAdapter, sql } from 'kysely'

const tables = ['auth_session', 'auth_account', 'auth_verification', 'auth_api_key']

export async function up(db: Kysely<any>): Promise<void> {
  if (!(db.getExecutor().adapter instanceof PostgresAdapter)) return
  for (const table of tables) {
    // (the columns are text, gen_random_uuid() is built into postgres 13+)
    await sql`alter table ${sql.table(table)} alter column "id" set default (gen_random_uuid())::text`.execute(db)
  }
}

export async function down(db: Kysely<any>): Promise<void> {
  if (!(db.getExecutor().adapter instanceof PostgresAdapter)) return
  for (const table of tables) {
    await sql`alter table ${sql.table(table)} alter column "id" drop default`.execute(db)
  }
}
