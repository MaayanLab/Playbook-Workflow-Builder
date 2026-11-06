import type { Adapter, AdapterAccount, AdapterSession, AdapterUser, VerificationToken } from 'next-auth/adapters'
import db from '@/app/db'
import { uuid, type Account, type Session, type User, type VerificationToken as VerificationTokenRecord } from '@/db'

const toDate = (value: string | null) => value !== null ? new Date(value) : null

function toAdapterUser(user: User): AdapterUser {
  return {
    ...user,
    name: user.name ?? undefined,
    image: user.image ?? undefined,
    emailVerified: toDate(user.emailVerified),
  }
}

function toAdapterSession(session: Session): AdapterSession {
  return { ...session, expires: toDate(session.expires) ?? new Date(0) }
}

function toVerificationToken(token: VerificationTokenRecord): VerificationToken {
  return {
    identifier: token.identifier ?? '',
    token: token.token ?? '',
    expires: toDate(token.expires) ?? new Date(0),
  }
}

function toAccountRecord(account: AdapterAccount) {
  return {
    userId: account.userId,
    type: account.type,
    provider: account.provider,
    providerAccountId: account.providerAccountId,
    refresh_token: account.refresh_token ?? null,
    access_token: account.access_token ?? null,
    expires_at: account.expires_at ?? null,
    token_type: account.token_type ?? null,
    scope: account.scope ?? null,
    id_token: account.id_token ?? null,
    session_state: account.session_state ?? null,
    oauth_token_secret: (account.oauth_token_secret as string | undefined) ?? null,
    oauth_token: (account.oauth_token as string | undefined) ?? null,
  } satisfies Omit<Account, 'id'>
}

export default function ORMAdapter(): Adapter {
  return {
    createUser: async (data) => {
      const user = await db.insertInto('user')
        .values({
          id: uuid(),
          name: data.name ?? null,
          email: data.email ?? '',
          emailVerified: data.emailVerified?.toISOString() ?? null,
          image: data.image ?? null,
        })
        .returningAll()
        .executeTakeFirstOrThrow()
      return toAdapterUser(user)
    },
    getUser: async (id) => {
      const user = await db.selectFrom('user').selectAll().where('id', '=', id).executeTakeFirst()
      return user ? toAdapterUser(user) : null
    },
    getUserByEmail: async (email) => {
      if (!email) throw new Error('Cannot search by email without one')
      const user = await db.selectFrom('user').selectAll().where('email', '=', email).executeTakeFirst()
      return user ? toAdapterUser(user) : null
    },
    getUserByAccount: async ({ provider, providerAccountId }) => {
      const user = await db.selectFrom('account')
        .innerJoin('user', 'user.id', 'account.userId')
        .selectAll('user')
        .where('account.provider', '=', provider)
        .where('account.providerAccountId', '=', providerAccountId)
        .executeTakeFirst()
      return user ? toAdapterUser(user) : null
    },
    updateUser: async ({ id, ...data_ }) => {
      const data: Partial<Omit<User, 'id' | 'created'>> = {}
      if (data_.name !== undefined) data.name = data_.name ?? ''
      if (data_.email !== undefined) data.email = data_.email
      if (data_.emailVerified !== undefined) data.emailVerified = data_.emailVerified?.toISOString() ?? null
      if (data_.image !== undefined) data.image = data_.image ?? ''
      const query = Object.keys(data).length > 0
        ? db.updateTable('user').set(data).where('id', '=', id).returningAll()
        : db.selectFrom('user').selectAll().where('id', '=', id)
      const user = await query.executeTakeFirst()
      if (!user) throw new Error('User not found')
      return toAdapterUser(user)
    },
    deleteUser: async (id) => {
      await db.deleteFrom('user').where('id', '=', id).execute()
    },
    linkAccount: async (account) => {
      await db.insertInto('account').values({ id: uuid(), ...toAccountRecord(account) }).execute()
    },
    unlinkAccount: async ({ provider, providerAccountId }) => {
      await db.deleteFrom('account')
        .where('provider', '=', provider)
        .where('providerAccountId', '=', providerAccountId)
        .execute()
    },
    getSessionAndUser: async (sessionToken) => {
      const session = await db.selectFrom('session').selectAll().where('sessionToken', '=', sessionToken).executeTakeFirst()
      if (!session) return null
      const user = await db.selectFrom('user').selectAll().where('id', '=', session.userId).executeTakeFirst()
      if (!user) return null
      return { session: toAdapterSession(session), user: toAdapterUser(user) }
    },
    createSession: async (data) => {
      const session = await db.insertInto('session')
        .values({ id: uuid(), sessionToken: data.sessionToken, userId: data.userId, expires: data.expires.toISOString() })
        .returningAll()
        .executeTakeFirstOrThrow()
      return toAdapterSession(session)
    },
    updateSession: async ({ sessionToken, ...data }) => {
      const session = await db.updateTable('session')
        .set({
          ...(data.userId !== undefined ? { userId: data.userId } : {}),
          ...(data.expires !== undefined ? { expires: data.expires.toISOString() } : {}),
          sessionToken,
        })
        .where('sessionToken', '=', sessionToken)
        .returningAll()
        .executeTakeFirst()
      return session ? toAdapterSession(session) : null
    },
    deleteSession: async (sessionToken) => {
      await db.deleteFrom('session').where('sessionToken', '=', sessionToken).execute()
    },
    createVerificationToken: async (data) => {
      const token = await db.insertInto('verification_token')
        .values({ id: uuid(), identifier: data.identifier, token: data.token, expires: data.expires.toISOString() })
        .returningAll()
        .executeTakeFirstOrThrow()
      return toVerificationToken(token)
    },
    useVerificationToken: async ({ identifier, token }) => {
      const deleted = await db.deleteFrom('verification_token')
        .where('identifier', '=', identifier)
        .where('token', '=', token)
        .returningAll()
        .executeTakeFirst()
      return deleted ? toVerificationToken(deleted) : null
    },
  }
}
