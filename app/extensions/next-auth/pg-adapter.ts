import * as pg from 'pg'
import * as db from '@/db/accounts'
import type { Adapter, AdapterSession, AdapterUser, VerificationToken } from 'next-auth/adapters'
import { ClientHelper } from '@/utils/orm'

export default function PgAdapter(connectionString: string): Adapter {
  const pool = new pg.Pool({ connectionString })
  const poolHelper = new ClientHelper(pool)
  return {
    createUser: async (data) => {
      return (await poolHelper.insertReturning(db.user, data)).one()
    },
    getUser: async (id) => {
      return (await poolHelper.selectWhere(db.user, { id })).oneOrNone()
    },
    getUserByEmail: async (email) => {
      const results = await poolHelper.selectWhere(db.user, { email })
      return results.oneOrNone()
    },
    getUserByAccount: async (provider_providerAccountId) => {
      const results = await poolHelper.query(subst => `
        select user.*
        from account
        inner join user on account.userId = user.id
        where account.provider = ${subst(provider_providerAccountId.provider)}
        and account.providerAccountId = ${subst(provider_providerAccountId.providerAccountId)};
      `, db.user.codec)
      return results.oneOrNone()
    },
    updateUser: async (data) => {
      if (data.id === undefined) throw new Error('NotImplementedError')
      return (await poolHelper.updateReturning(db.user, data)).oneOrNone()
    },
    deleteUser: async (id) => {
      await poolHelper.deleteWhere(db.user, { id })
    },
    linkAccount: async (data) => {
      await poolHelper.insert(db.account, data)
    },
    unlinkAccount: async (provider_providerAccountId) => {
      await poolHelper.deleteWhere(db.account, provider_providerAccountId)
    },
    async getSessionAndUser(sessionToken) {
      const results = await poolHelper.query(subst => `
        select
          "session"."id" as "session_id",
          "session"."expires" as "session_expires",
          "session"."sessionToken" as "session_sessionToken",
          "session"."userId" as "session_userId",
          "user"."id" as "user_id",
          "user"."name" as "user_name",
          "user"."email" as "user_email",
          "user"."emailVerified" as "user_emailVerified",
          "user"."image" as "user_image"
        from "session"
        inner join "user" on "session"."userId" = "user"."id"
        where "session"."sessionToken" = ${subst(sessionToken)};
      `)
      const result = results.one()
      if (result === null) return null
      else {
        const session: AdapterSession = db.session.codec.decode({
          id: result.session_id,
          expires: result.session_expires,
          sessionToken: result.session_sessionToken,
          userId: result.session_userId,
        })
        const user: AdapterUser = db.user.codec.decode({
          id: result.user_id,
          name: result.user_name,
          email: result.user_email,
          emailVerified: result.user_emailVerified,
          image: result.user_image,
          affiliation: '',
        })
        return { session, user }
      }
    },
    createSession: async (data) => {
      await poolHelper.insert(db.session, data)
      return data
    },
    updateSession: async (data) => {
      if (data.id === undefined) throw new Error('NotImplementedError')
      return (await poolHelper.updateReturning(db.session, data)).oneOrNone()
    },
    deleteSession: async (sessionToken) => {
      await poolHelper.deleteWhere(db.session, { sessionToken })
    },
    createVerificationToken: async (data) => {
      await poolHelper.insert(db.verification_token, data)
      return data
    },
    useVerificationToken: async (identifier_token) => {
      return (await poolHelper.deleteWhereReturning(db.verification_token, identifier_token)).oneOrNone()
    },
  }
}
