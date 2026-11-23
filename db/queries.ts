/**
 * Reusable queries
 */
import type { Insertable } from 'kysely'
import { type DB, uuid } from '.'
import type { Database } from './schema'
import { tsvector } from '../utils/tsvector'
import { type UserPlaybookSearchable, userPlaybookIndexRows } from './user_playbook_index'
export { splitTypes } from './user_playbook_index'

/** user_upload joined with its upload */
export function selectUserUploadComplete(db: DB) {
  return db.selectFrom('user_upload')
    .innerJoin('upload', 'upload.id', 'user_upload.upload')
    .select([
      'user_upload.id',
      'user_upload.user',
      'upload.url',
      'upload.sha256',
      'upload.size',
      'user_upload.filename',
      'user_upload.created',
    ])
}

/** Record an upload (url & checksum) for a user (or anonymously), returning the user_upload id */
export async function insertUserUpload(db: DB, file: { url: string, sha256: string, size: number, filename: string }, user: string | null) {
  return await db.transaction().execute(async (trx) => {
    const upload = await trx.insertInto('upload')
      .values({ id: uuid(), url: file.url, sha256: file.sha256, size: file.size })
      .returning('id')
      .executeTakeFirstOrThrow()
    return await trx.insertInto('user_upload')
      .values({ id: uuid(), user, upload: upload.id, filename: file.filename })
      .returning('id')
      .executeTakeFirstOrThrow()
  })
}

/* ---------------------------------------------------------------- user_playbook search */

async function insertChunked<T extends 'user_playbook_trigram' | 'user_playbook_type'>(db: DB, table: T, rows: Insertable<Database[T]>[]) {
  // stay well under sqlite's bound parameter limit
  for (let i = 0; i < rows.length; i += 500) {
    await db.insertInto(table).values(rows.slice(i, i + 500) as any).execute()
  }
}

/** (Re)build the search index for a user_playbook, call this whenever its text or types change */
export async function indexUserPlaybook(db: DB, id: string, record: UserPlaybookSearchable) {
  const { trigrams, types } = userPlaybookIndexRows(id, record)
  await db.deleteFrom('user_playbook_trigram').where('user_playbook', '=', id).execute()
  await db.deleteFrom('user_playbook_type').where('user_playbook', '=', id).execute()
  await insertChunked(db, 'user_playbook_trigram', trigrams)
  await insertChunked(db, 'user_playbook_type', types)
}

/**
 * Filter (by input/output types) & rank (by trigram overlap with the search, then clicks)
 *  user_playbooks in the db. The score is the number of the search's trigrams the playbook has.
 */
export function searchUserPlaybooks(db: DB, props: { search?: string, inputs?: string[], outputs?: string[] }) {
  let query = db.selectFrom('user_playbook').selectAll('user_playbook')
  for (const [kind, specs] of [['input', props.inputs], ['output', props.outputs]] as const) {
    for (const spec of specs ?? []) {
      query = query.where(eb => eb.exists(
        eb.selectFrom('user_playbook_type')
          .select('user_playbook_type.user_playbook')
          .whereRef('user_playbook_type.user_playbook', '=', 'user_playbook.id')
          .where('user_playbook_type.kind', '=', kind)
          .where('user_playbook_type.spec', '=', spec)
      ))
    }
  }
  const trigrams = props.search ? [...tsvector(props.search)] : []
  if (trigrams.length === 0) return query.orderBy('user_playbook.clicks', 'desc')
  return query
    .innerJoin(
      eb => eb.selectFrom('user_playbook_trigram')
        .select(['user_playbook_trigram.user_playbook', eb.fn.countAll().as('score')])
        .where('user_playbook_trigram.trigram', 'in', trigrams)
        .groupBy('user_playbook_trigram.user_playbook')
        .as('match'),
      join => join.onRef('match.user_playbook', '=', 'user_playbook.id'),
    )
    .orderBy('match.score', 'desc')
    .orderBy('user_playbook.clicks', 'desc')
}
