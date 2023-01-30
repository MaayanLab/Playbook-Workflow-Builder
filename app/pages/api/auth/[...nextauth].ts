import NextAuth, { NextAuthOptions, Session } from 'next-auth'
import type { Provider } from "next-auth/providers"
import EmailProvider from 'next-auth/providers/email'
import GoogleProvider from "next-auth/providers/google"
import MemoryAdapter from '@/app/extensions/next-auth/memory-adapter'
import PgAdapter from '@/app/extensions/next-auth/pg-adapter'

export type AuthOptions = typeof authOptions
export type SessionWithId = Session & { user: { id: string } }

export const authOptions: NextAuthOptions = {
  adapter: process.env.DATABASE_URL ? PgAdapter(process.env.DATABASE_URL) : MemoryAdapter(),
  providers: ([
    process.env.NEXTAUTH_GOOGLE ? GoogleProvider(JSON.parse(process.env.NEXTAUTH_GOOGLE)) : undefined,
    process.env.EMAIL_SERVER && process.env.EMAIL_FROM ? EmailProvider({
      server: process.env.EMAIL_SERVER,
      from: process.env.EMAIL_FROM,
    }) : undefined,
  ] as Array<Provider | undefined>).filter((v): v is Provider => v !== undefined),
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
