/**
 * Exports (LaTeX bundles, report PDFs, ...) of a playbook which are produced in the background
 *  by workers, rather than in the request (which can exceed the ingress timeout). Exports belong
 *  to the (signed in) user who requested them, the result is stored as one of their uploads.
 */
import { type Kysely, PostgresAdapter, sql } from 'kysely'

export async function up(db: Kysely<any>): Promise<void> {
  const pg = db.getExecutor().adapter instanceof PostgresAdapter
  const uuid = pg ? 'uuid' : 'text'
  const timestamp = pg ? 'timestamp' : 'text'
  const now = pg ? sql`now()` : sql`(strftime('%Y-%m-%dT%H:%M:%fZ', 'now'))`
  await db.schema.createTable('fpl_export')
    .addColumn('id', uuid, col => col.primaryKey())
    .addColumn('fpl', uuid, col => col.notNull().references('fpl.id').onDelete('cascade'))
    // who requested the export, they own it & the file it produces
    .addColumn('user', uuid, col => col.notNull().references('user.id').onDelete('cascade'))
    .addColumn('format', 'text', col => col.notNull())
    // format specific parameters (JSON), e.g. report metadata & author
    .addColumn('params', 'text', col => col.notNull())
    // queued | running | complete | error
    .addColumn('status', 'text', col => col.notNull())
    .addColumn('error', 'text')
    // the produced file, once complete
    .addColumn('upload', uuid, col => col.references('user_upload.id').onDelete('set null'))
    .addColumn('created', timestamp, col => col.notNull().defaultTo(now))
    .addColumn('updated', timestamp, col => col.notNull().defaultTo(now))
    .execute()
  await db.schema.createIndex('fpl_export_user_idx').on('fpl_export').column('user').execute()
}

export async function down(db: Kysely<any>): Promise<void> {
  await db.schema.dropTable('fpl_export').execute()
}
