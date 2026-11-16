/**
 * thread.openai_thread was `not null` in databases provisioned with dbmate,
 *  but threads are created without one. (sqlite databases already have it nullable)
 */
import { type Kysely, PostgresAdapter, sql } from 'kysely'

export async function up(db: Kysely<any>): Promise<void> {
  if (db.getExecutor().adapter instanceof PostgresAdapter) {
    await sql`alter table "thread" alter column "openai_thread" drop not null`.execute(db)
  }
}

export async function down(db: Kysely<any>): Promise<void> {
  // irreversible: threads without an openai_thread would violate the constraint
}
