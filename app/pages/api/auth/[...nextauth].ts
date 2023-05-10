import NextAuth, { ISODateString, NextAuthOptions, Session } from 'next-auth'
import type { Provider } from "next-auth/providers"
import db from '@/app/db'
import CredentialsProvider from 'next-auth/providers/credentials'
import EmailProvider from 'next-auth/providers/email'
import GoogleProvider from "next-auth/providers/google"
import ORMAdapter from '@/app/extensions/next-auth/orm-adapter'
import ORCIDProvider from '@/app/extensions/next-auth/orcid-provider'

export type AuthOptions = typeof authOptions
export type SessionWithId = {
  user: Session['user'] & { id: string },
  expires: ISODateString,
}

export const authOptions: NextAuthOptions = {
  adapter: ORMAdapter(),
  providers: ([
    process.env.NODE_ENV === 'development' ? CredentialsProvider({
      id: 'dev',
      name: 'Development User',
      credentials: {},
      async authorize(credentials, req) {
        return await db.objects.user.upsert({
          where: {
            id: '00000000-0000-0000-0000-000000000000',
          },
          create: {
            id: '00000000-0000-0000-0000-000000000000',
            name: 'Development User',
            email: '',
          },
        })
      }
    }) : undefined,
    process.env.NEXTAUTH_GOOGLE ? GoogleProvider(JSON.parse(process.env.NEXTAUTH_GOOGLE)) : undefined,
    process.env.NEXTAUTH_ORCID ? ORCIDProvider(JSON.parse(process.env.NEXTAUTH_ORCID), db) : undefined,
    process.env.EMAIL_SERVER && process.env.EMAIL_FROM ? EmailProvider({
      server: process.env.EMAIL_SERVER,
      from: process.env.EMAIL_FROM,
    }) : undefined,
  ] as Array<Provider | undefined>).filter((v): v is Provider => v !== undefined),
  callbacks: {
    session: async ({ session, token }) => {
      Object.assign(session, {
        user: await db.objects.user.findUnique({
          where: {
            id: token.sub,
          }
        })
      })
      return session as SessionWithId
    }
  },
  session: {
    strategy: 'jwt',
  },
}

export default NextAuth(authOptions)
