/**
 * Preserve the typed Responses API item behind each displayed chat message.
 * `content` and `role` remain the UI projection and support legacy threads.
 */
import type { Kysely } from 'kysely'

export async function up(db: Kysely<any>): Promise<void> {
  await db.schema.alterTable('thread_message')
    .addColumn('response_item', 'text')
    .execute()
}

export async function down(db: Kysely<any>): Promise<void> {
  await db.schema.alterTable('thread_message').dropColumn('response_item').execute()
}
