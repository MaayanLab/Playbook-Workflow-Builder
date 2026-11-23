/**
 * Search indices for user_playbooks so that community playbooks can be filtered (by input/output
 *  data types) & ranked (by trigram overlap with the search, see utils/tsvector) in the db, rather
 *  than by loading every public playbook. Kept up to date by indexUserPlaybook (db/queries.ts).
 */
import { type Kysely, PostgresAdapter } from 'kysely'
import { userPlaybookIndexRows } from '../user_playbook_index'

export async function up(db: Kysely<any>): Promise<void> {
  const pg = db.getExecutor().adapter instanceof PostgresAdapter
  const uuid = pg ? 'uuid' : 'text'
  await db.schema.createTable('user_playbook_trigram')
    .addColumn('user_playbook', uuid, col => col.notNull().references('user_playbook.id').onDelete('cascade'))
    .addColumn('trigram', 'text', col => col.notNull())
    .addPrimaryKeyConstraint('user_playbook_trigram_pkey', ['user_playbook', 'trigram'])
    .execute()
  await db.schema.createIndex('user_playbook_trigram_trigram_idx').on('user_playbook_trigram').columns(['trigram', 'user_playbook']).execute()
  await db.schema.createTable('user_playbook_type')
    .addColumn('user_playbook', uuid, col => col.notNull().references('user_playbook.id').onDelete('cascade'))
    .addColumn('kind', 'text', col => col.notNull())
    .addColumn('spec', 'text', col => col.notNull())
    .addPrimaryKeyConstraint('user_playbook_type_pkey', ['user_playbook', 'kind', 'spec'])
    .execute()
  await db.schema.createIndex('user_playbook_type_spec_idx').on('user_playbook_type').columns(['kind', 'spec', 'user_playbook']).execute()

  // backfill
  let after: string | undefined
  while (true) {
    let query = db.selectFrom('user_playbook').select(['id', 'title', 'description', 'inputs', 'outputs']).orderBy('id').limit(500)
    if (after !== undefined) query = query.where('id', '>', after)
    const playbooks = await query.execute()
    if (playbooks.length === 0) break
    after = playbooks[playbooks.length - 1].id
    const trigrams = [], types = []
    for (const playbook of playbooks) {
      const rows = userPlaybookIndexRows(playbook.id, playbook)
      trigrams.push(...rows.trigrams)
      types.push(...rows.types)
    }
    for (let i = 0; i < trigrams.length; i += 500) await db.insertInto('user_playbook_trigram').values(trigrams.slice(i, i + 500)).execute()
    for (let i = 0; i < types.length; i += 500) await db.insertInto('user_playbook_type').values(types.slice(i, i + 500)).execute()
  }
}

export async function down(db: Kysely<any>): Promise<void> {
  await db.schema.dropTable('user_playbook_type').execute()
  await db.schema.dropTable('user_playbook_trigram').execute()
}
