import NextAuth, { NextAuthOptions, Session } from 'next-auth'
import EmailProvider from 'next-auth/providers/email'
import MemoryAdapter from '@/app/extensions/next-auth/memory-adapter'
import PgAdapter from '@/app/extensions/next-auth/pg-adapter'

export type AuthOptions = typeof authOptions
export type SessionWithId = Session & { user: { id: string } }

export const authOptions: NextAuthOptions = {
  adapter: process.env.DATABASE_URL ? PgAdapter(process.env.DATABASE_URL) : MemoryAdapter(),
  providers: [
    EmailProvider({
      server: process.env.EMAIL_SERVER,
      from: process.env.EMAIL_FROM,
    })
  ],
  callbacks: {
    session: async ({ session, token }) => {
      if (session?.user) Object.assign(session.user, { id: token.sub })
      return session as SessionWithId
    }
  },
  session: {
    strategy: 'jwt',
  },
}

export default NextAuth(authOptions)
