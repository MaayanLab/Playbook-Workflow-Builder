import NextAuth, { NextAuthOptions, Session } from 'next-auth'
import type { Provider } from "next-auth/providers"
import EmailProvider from 'next-auth/providers/email'
import GoogleProvider from "next-auth/providers/google"
import ORMAdapter from '@/app/extensions/next-auth/orm-adapter'
import ORCIDProvider from '@/app/extensions/next-auth/orcid-provider'

export type AuthOptions = typeof authOptions
export type SessionWithId = Session & { user: { id: string } }

export const authOptions: NextAuthOptions = {
  adapter: ORMAdapter(),
  providers: ([
    process.env.NEXTAUTH_GOOGLE ? GoogleProvider(JSON.parse(process.env.NEXTAUTH_GOOGLE)) : undefined,
    process.env.NEXTAUTH_ORCID ? ORCIDProvider(JSON.parse(process.env.NEXTAUTH_ORCID)) : undefined,
    process.env.EMAIL_SERVER && process.env.EMAIL_FROM ? EmailProvider({
      server: process.env.EMAIL_SERVER,
      from: process.env.EMAIL_FROM,
    }) : undefined,
  ] as Array<Provider | undefined>).filter((v): v is Provider => v !== undefined),
  callbacks: {
    session: async ({ session, user }) => {
      if (session?.user && user) Object.assign(session.user, { id: user.id })
      return session as SessionWithId
    }
  },
  session: {
    strategy: 'database',
  },
}

export default NextAuth(authOptions)
