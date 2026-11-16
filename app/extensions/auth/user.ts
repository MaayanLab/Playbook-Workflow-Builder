import db from '@/app/db'

export async function getUser(id: string) {
  return await db.selectFrom('user').selectAll().where('id', '=', id).executeTakeFirst()
}

export async function getUserOrcidAccount(userId: string) {
  const account = await db.selectFrom('auth_account').selectAll()
    .where('userId', '=', userId)
    .where('providerId', '=', 'orcid')
    .executeTakeFirst()
  return account ? {
    ...account,
    providerAccountId: account.accountId,
    provider: account.providerId,
    access_token: account.accessToken,
    refresh_token: account.refreshToken,
    id_token: account.idToken,
    expires_at: account.accessTokenExpiresAt ? new Date(account.accessTokenExpiresAt).getTime() / 1000 : null,
  } : undefined
}

