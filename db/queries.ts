/**
 * Reusable queries
 */
import { type DB, uuid } from '.'

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
