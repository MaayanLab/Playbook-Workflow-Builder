import * as pg from 'pg'
import * as db from '@/db/accounts'
import type { Adapter, AdapterSession, AdapterUser, VerificationToken } from 'next-auth/adapters'

export default function PgAdapter(connectionString: string): Adapter {
  const pool = new pg.Pool({ connectionString })
  return {
    createUser: async (data) => {
      const results = await pool.query(`
        insert into user ("name", "email", "emailVerified", "image")
        values ($1, $2, $3, $4) returning id
      `, [data.name, data.email, data.emailVerified, data.image])
      const [[id]] = results.rows
      return { ...data, id }
    },
    getUser: async (id) => {
      try {
        const results = await pool.query(
          'select * from user where "id" = $1',
          [id])
        if (results.rowCount === 0) return null
        else {
          const user: AdapterUser = db.user.codec.decode(results.rows[0])
          return user
        }
      } catch (e) {
        return null
      }
    },
    getUserByEmail: async (email) => {
      try {
        const results = await pool.query(`
          select * from user where "email" = $1;
        `, [email])
        if (results.rowCount === 0) return null
        else {
          const user: AdapterUser = db.user.codec.decode(results.rows[0])
          return user
        }
      } catch (e) {
        return null
      }
    },
    getUserByAccount: async (provider_providerAccountId) => {
      const results = await pool.query(`
        select user.*
        from account
        inner join user on account.userId = user.id
        where account.provider = $1 and account.providerAccountId = $2;
        `, [provider_providerAccountId.provider, provider_providerAccountId.providerAccountId])
      if (results.rowCount === 0) return null
      else {
        const user: AdapterUser = db.user.codec.decode(results.rows[0])
        return user
      }
    },
    updateUser: async (data) => {
      const results = await pool.query(`
        update user
        set "name" = coalesce($2, "name"),
            "email" = coalesce($3, "email"),
            "emailVerified" = coalesce($4, "emailVerified"),
            "image" = coalesce($5, "image")
        where "id" = $1
        returning *;
      `, [data.id, data.name, data.email, data.emailVerified, data.image])
      if (results.rowCount === 0) throw new Error('User not found')
      else {
        const user: AdapterUser = db.user.codec.decode(results.rows[0])
        return user
      }
    },
    deleteUser: async (id) => {
      await pool.query(`
        delete from user
        where "id" = $1;
      `, [id])
    },
    linkAccount: async (data) => {
      await pool.query(`
        insert into account ("userId", "type", "provider", "providerAccountId", "refresh_token", "access_token", "expires_at", "token_type", "scope", "id_token", "session_state", "oauth_token_secret", "oauth_token")
        values ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13);
      `, [data.userId, data.type, data.provider, data.providerAccountId, data.refresh_token, data.access_token, data.expires_at, data.token_type, data.scope, data.id_token, data.session_state, data.oauth_token_secret, data.oauth_token])
    },
    unlinkAccount: async (provider_providerAccountId) => {
      await pool.query(`
        delete from account
        where "provider" = $1 and providerAccountId = $2;
      `, [provider_providerAccountId.provider, provider_providerAccountId.providerAccountId])
    },
    async getSessionAndUser(sessionToken) {
      const results = await pool.query(`
        select
          row_to_jsonb(session.*) as session,
          row_to_jsonb(user.*) as user
        from session
        inner join user on session.userId = user.id
        where session.sessionToken = $1;
        `, [sessionToken])
      if (results.rowCount === 0) return null
      else {
        const session: AdapterSession = db.session.codec.decode(results.rows[0].session)
        const user: AdapterUser = db.user.codec.decode(results.rows[0].user)
        return { session, user }
      }
    },
    createSession: async (data) => {
      await pool.query(`
        insert into session ("expires", "sessionToken", "userId")
        values ($1, $2, $3);
      `, [data.expires, data.sessionToken, data.userId])
      return data
    },
    updateSession: async (data) => {
      const results = await pool.query(`
        update session
        set "userId" = coalesce($2, "userId"),
            "expires" = coalesce($3, "expires")
        where "sessionToken" = $1
        returning *;
      `, [data.sessionToken, data.userId, data.expires])
      if (results.rowCount === 0) return null
      else {
        const session: AdapterSession = db.session.codec.decode(results.rows[0])
        return session
      }
    },
    deleteSession: async (sessionToken) => {
      await pool.query(`
        delete from session
        where "sessionToken" = $1;
      `, [sessionToken])
    },
    createVerificationToken: async (data) => {
      await pool.query(`
        insert into verification_token ("identifier", "token", "expires")
        values ($1, $2, $3);
      `, [data.identifier, data.token, data.expires])
      return data
    },
    useVerificationToken: async (identifier_token) => {
      const results = await pool.query(`
        delete from verification_token
        where "identifier" = $1 and "token" = $2
        returning *;
      `, [identifier_token.identifier, identifier_token.token])
      if (results.rowCount === 0) return null
      else {
        const verificationToken: VerificationToken = db.verification_token.codec.decode(results.rows[0])
        return verificationToken
      }
    },
  }
}
