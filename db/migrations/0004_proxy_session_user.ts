/**
 * Record which user launched a CAVATICA session, so that only they can use it
 *  and its worker can authenticate with their CAVATICA API key.
 *  (sessions are ephemeral, existing ones are left without a user)
 */
import { type Kysely, PostgresAdapter } from 'kysely'

export async function up(db: Kysely<any>): Promise<void> {
  const pg = db.getExecutor().adapter instanceof PostgresAdapter
  await db.schema.alterTable('proxy_session')
    .addColumn('user', pg ? 'uuid' : 'text', col => col.references('user.id').onDelete('cascade'))
    .execute()
}

export async function down(db: Kysely<any>): Promise<void> {
  await db.schema.alterTable('proxy_session').dropColumn('user').execute()
}
