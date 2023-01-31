import type { Adapter } from 'next-auth/adapters'
import db from '@/app/db'

export default function ORMAdapter(): Adapter {
  return {
    createUser: async (data) => {
      return await db.objects.user.create({ data })
    },
    getUser: async (id) => {
      return await db.objects.user.findUnique({ where: { id } })
    },
    getUserByEmail: async (email) => {
      return await db.objects.user.findUnique({ where: { email } })
    },
    getUserByAccount: async (provider_providerAccountId) => {
      // TODO: return await db.objects.user.findUnique({ where: { accounts: provider_providerAccountId } })
      // @ts-ignore
      const account = await db.objects.account.findUnique({ where: provider_providerAccountId })
      const user = (account !== null) ? await db.objects.user.findUnique({ where: { id: account.userId } }) : null
      return user
    },
    updateUser: async (data) => {
      const user = await db.objects.user.update({ where: { id: data.id }, data })
      if (user === null) throw new Error('User not found')
      return user
    },
    deleteUser: async (id) => {
      await db.objects.user.delete({ where: { id } })
    },
    linkAccount: async (data) => {
      await db.objects.account.create({ data })
    },
    unlinkAccount: async (provider_providerAccountId) => {
      await db.objects.account.delete({ where: provider_providerAccountId })
    },
    getSessionAndUser: async (sessionToken) => {
      // const result = await db.objects.session.findUnique({ where: { sessionToken }, include: { user: true } })
      const session = await db.objects.session.findUnique({ where: { sessionToken } })
      const user = (session !== null) ? await db.objects.user.findUnique({ where: { id: session.userId } }) : null
      return session !== null && user !== null ? { session, user } : null
    },
    createSession: async (data) => {
      return await db.objects.session.create({ data })
    },
    updateSession: async (data) => {
      return await db.objects.session.update({ where: { sessionToken: data.sessionToken }, data })
    },
    deleteSession: async (sessionToken) => {
      await db.objects.session.delete({ where: { sessionToken } })
    },
    createVerificationToken: async (data) => {
      return await db.objects.verification_token.create({ data })
    },
    useVerificationToken: async (identifier_token) => {
      return await db.objects.verification_token.delete({ where: identifier_token })
    },
  }
}