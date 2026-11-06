import db from '@/app/db'

export async function getUser(id: string) {
  return await db.selectFrom('user').selectAll().where('id', '=', id).executeTakeFirst()
}

export async function getUserOrcidAccount(userId: string) {
  return await db.selectFrom('account').selectAll()
    .where('userId', '=', userId)
    .where('provider', '=', 'orcid')
    .executeTakeFirst()
}
