import { LevelUp } from 'levelup'
import type { Adapter, AdapterAccount, AdapterSession, AdapterUser, VerificationToken } from 'next-auth/adapters'
import { z } from 'zod'
import sub from 'subleveldown'

const DateC = z.union([z.date(), z.string()]).transform(date => new Date(date))
const AdapterUserC = z.object({
  id: z.string(),
  email: z.string(),
  emailVerified: z.union([DateC, z.null()]),
})
const AdapterSessionC = z.object({
  sessionToken: z.string(),
  userId: z.string(),
  expires: DateC,
})

export default function LevelAdapter(db: LevelUp): Adapter {
  const dbNextAuth: LevelUp = sub(db, 'nextauth')
  const dbUser: LevelUp = sub(dbNextAuth, 'user')
  const dbAccount: LevelUp = sub(dbNextAuth, 'account')
  const dbSession: LevelUp = sub(dbNextAuth, 'session')
  const dbVerificationToken: LevelUp = sub(dbNextAuth, 'verificationToken')
  return {
    createUser: async (data) => {
      const user = { ...data, id: data.email }
      await dbUser.put(user.id, JSON.stringify(user))
      return user
    },
    getUser: async (id) => {
      try {
        const user: AdapterUser = AdapterUserC.parse(JSON.parse(await dbUser.get(`${id}`)))
        return user
      } catch (e) {
        return null
      }
    },
    getUserByEmail: async (email) => {
      try {
        const user: AdapterUser = AdapterUserC.parse(JSON.parse(await dbUser.get(`${email}`)))
        return user
      } catch (e) {
        return null
      }
    },
    getUserByAccount: async (provider_providerAccountId) => {
      let account: AdapterAccount | null, user: AdapterUser | null
      try {
        account = JSON.parse(await dbAccount.get(`${provider_providerAccountId}`))
      } catch (e) {
        account = null
      }
      try {
        user = account ? AdapterUserC.parse(JSON.parse(await dbUser.get(`${account.userId}`))) : null
      } catch (e) {
        user = null
      }
      return user
    },
    updateUser: async (data) => {
      const user: AdapterUser = AdapterUserC.parse(JSON.parse(await dbUser.get(`${data.id}`)))
      Object.assign(user, data)
      await dbUser.put(`${data.id}`, JSON.stringify(user))
      return user
    },
    deleteUser: async (id) => {
      await dbUser.del(`${id}`)
    },
    linkAccount: async (data) => {
      await dbAccount.put(`${data.provider}_${data.providerAccountId}`, JSON.stringify(data))
    },
    unlinkAccount: async (provider_providerAccountId) => {
      await dbAccount.del(`${provider_providerAccountId}`)
    },
    async getSessionAndUser(sessionToken) {
      let session: AdapterSession | null, user: AdapterUser | null
      try {
        session = AdapterSessionC.parse(JSON.parse(await dbSession.get(`${sessionToken}`)))
      } catch (e) {
        session = null
      }
      try {
        user = session ? AdapterUserC.parse(JSON.parse(await dbUser.get(`${session.userId}`))) : null
      } catch (e) {
        user = null
      }
      return user && session ? { user, session } : null
    },
    createSession: async (data) => {
      await dbSession.put(`${data.sessionToken}`, JSON.stringify(data))
      return data
    },
    updateSession: async (data) => {
      let session: AdapterSession = AdapterSessionC.parse(JSON.parse(await dbSession.get(`${data.sessionToken}`)))
      Object.assign(session, data)
      await dbSession.put(`${data.sessionToken}`, JSON.stringify(session))
      return session
    },
    deleteSession: async (sessionToken) => {
      await dbSession.del(`${sessionToken}`)
    },
    createVerificationToken: async (data) => {
      await dbVerificationToken.put(`${data.identifier}_${data.token}`, JSON.stringify(data))
      return data
    },
    useVerificationToken: async (identifier_token) => {
      let verificationToken: VerificationToken | null
      try {
        verificationToken = JSON.parse(await dbVerificationToken.get(`${identifier_token.identifier}_${identifier_token.token}`))
        await dbVerificationToken.del(`${identifier_token.identifier}`)
      } catch (e) {
        verificationToken = null
      }
      return verificationToken
    },
  }
}