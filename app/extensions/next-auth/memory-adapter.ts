import type { Adapter, AdapterAccount, AdapterSession, AdapterUser, VerificationToken } from 'next-auth/adapters'

export default function MemoryAdapter(): Adapter {
  const users: Record<string, AdapterUser> = {}
  const accounts: Record<string, AdapterAccount> = {}
  const sessions: Record<string, AdapterSession> = {}
  const verificationTokens: Record<string, VerificationToken> = {}
  return {
    createUser: async (data) => {
      const user = { ...data, id: data.email }
      users[user.id] = user
      return user
    },
    getUser: async (id) => {
      return users[id] || null
    },
    getUserByEmail: async (email) => {
      return users[email] || null
    },
    getUserByAccount: async (provider_providerAccountId) => {
      const account = accounts[`${provider_providerAccountId.provider}_${provider_providerAccountId.providerAccountId}`] || null
      const user = users[`${account.userId}`] || null
      return user
    },
    updateUser: async (data) => {
      const id = data.id || data.email
      if (id === undefined) throw new Error('')
      const user = users[id]
      Object.assign(user, data)
      return user
    },
    deleteUser: async (id) => {
      delete users[id]
    },
    linkAccount: async (data) => {
      accounts[`${data.provider}_${data.providerAccountId}`] = data
    },
    unlinkAccount: async (provider_providerAccountId) => {
      delete accounts[`${provider_providerAccountId.provider}_${provider_providerAccountId.providerAccountId}`]
    },
    async getSessionAndUser(sessionToken) {
      const session = sessions[sessionToken] || null
      const user = session ? users[session.userId] || null : null
      return user && session ? { user, session } : null
    },
    createSession: async (data) => {
      sessions[data.sessionToken] = data
      return data
    },
    updateSession: async (data) => {
      const session = sessions[data.sessionToken]
      Object.assign(session, data)
      return session
    },
    deleteSession: async (sessionToken) => {
      delete sessions[sessionToken]
    },
    createVerificationToken: async (data) => {
      verificationTokens[`${data.identifier}_${data.token}`] = data
      return data
    },
    useVerificationToken: async (identifier_token) => {
      const verificationToken = verificationTokens[`${identifier_token.identifier}_${identifier_token.token}`] || null
      if (verificationToken !== null) {
        delete verificationTokens[`${identifier_token.identifier}_${identifier_token.token}`]
      }
      return verificationToken
    },
  }
}