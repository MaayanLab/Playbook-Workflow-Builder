/**
 * The DRS objects this instance serves:
 * - `<sha256>`: content addressed files (verified, i.e. stored as `storage://<sha256>`). Files are
 *    given these urls so that identical files are identical inputs, and the processes using them
 *    can be reused (the process graph is content addressed).
 * - `<user upload id>`: an upload (e.g. of a file which isn't stored here, but at a url).
 */
import db from '@/app/db'
import { selectUserUploadComplete } from '@/db/queries'

export const isContentId = (id: string) => /^[0-9a-f]{64}$/.test(id)

/** `drs://{this instance}/{id}` */
export const drsUrl = (id: string) => `${(process.env.PUBLIC_URL || '').replace(/^https?:/, 'drs:').replace(/\/$/, '')}/${id}`

export type DrsObject = { id: string, url: string, sha256: string, size: number, filename: string, created: string }

export async function resolveDrsObject(id: string): Promise<DrsObject | undefined> {
  if (isContentId(id)) {
    const upload = await db.selectFrom('upload')
      .select(['url', 'sha256', 'size', 'created'])
      .where('sha256', '=', id)
      .where('url', '=', `storage://${id}`)
      .orderBy('created')
      .executeTakeFirst()
    // (the content may have been uploaded under several names, by several users: it's named by its hash)
    return upload ? { id, ...upload, filename: id } : undefined
  }
  return await selectUserUploadComplete(db).where('user_upload.id', '=', id).executeTakeFirst()
}
