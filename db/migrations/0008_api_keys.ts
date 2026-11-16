/**
 * API keys (Better Auth's api-key plugin, see app/auth.ts), so that users can use the API
 *  (e.g. to upload files) as themselves. Only a hash of each key is stored.
 */
import { type Kysely, PostgresAdapter, sql } from 'kysely'

export async function up(db: Kysely<any>): Promise<void> {
  const pg = db.getExecutor().adapter instanceof PostgresAdapter
  const timestamp = pg ? 'timestamp' : 'text'
  const now = pg ? sql`now()` : sql`(strftime('%Y-%m-%dT%H:%M:%fZ', 'now'))`
  await db.schema.createTable('auth_api_key')
    .addColumn('id', 'text', col => col.primaryKey())
    .addColumn('configId', 'text', col => col.notNull().defaultTo('default'))
    .addColumn('name', 'text')
    // the first few characters of the key, to recognize it by
    .addColumn('start', 'text')
    .addColumn('prefix', 'text')
    // the (hashed) key
    .addColumn('key', 'text', col => col.notNull())
    .addColumn('referenceId', pg ? 'uuid' : 'text', col => col.notNull().references('user.id').onDelete('cascade'))
    .addColumn('refillInterval', 'bigint')
    .addColumn('refillAmount', 'bigint')
    .addColumn('lastRefillAt', timestamp)
    .addColumn('enabled', 'boolean', col => col.defaultTo(pg ? true : 1))
    .addColumn('rateLimitEnabled', 'boolean', col => col.defaultTo(pg ? true : 1))
    .addColumn('rateLimitTimeWindow', 'bigint')
    .addColumn('rateLimitMax', 'bigint')
    .addColumn('requestCount', 'bigint', col => col.defaultTo(0))
    .addColumn('remaining', 'bigint')
    .addColumn('lastRequest', timestamp)
    .addColumn('expiresAt', timestamp)
    .addColumn('createdAt', timestamp, col => col.notNull().defaultTo(now))
    .addColumn('updatedAt', timestamp, col => col.notNull().defaultTo(now))
    .addColumn('permissions', 'text')
    .addColumn('metadata', 'text')
    .execute()
  await db.schema.createIndex('auth_api_key_key_idx').on('auth_api_key').column('key').execute()
  await db.schema.createIndex('auth_api_key_reference_idx').on('auth_api_key').column('referenceId').execute()
}

export async function down(db: Kysely<any>): Promise<void> {
  await db.schema.dropTable('auth_api_key').execute()
}
