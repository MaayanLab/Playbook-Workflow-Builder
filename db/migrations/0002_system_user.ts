/**
 * Create the built-in user which owns application-created records (see @/db/constants)
 */
import type { Kysely } from 'kysely'
import { SYSTEM_USER_ID } from '../constants'

export async function up(db: Kysely<any>): Promise<void> {
  await db.insertInto('user')
    .values({ id: SYSTEM_USER_ID, name: 'Playbook Workflow Builder', email: '' })
    .onConflict(oc => oc.column('id').doNothing())
    .execute()
}

export async function down(db: Kysely<any>): Promise<void> {
  // the user may own records (and may predate this migration), so it is left in place
}
