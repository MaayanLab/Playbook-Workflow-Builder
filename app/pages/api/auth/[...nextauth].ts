import NextAuth from 'next-auth'
import type { Provider } from "next-auth/providers"
import EmailProvider from 'next-auth/providers/email'
import GoogleProvider from "next-auth/providers/google"
import MemoryAdapter from '@/app/extensions/next-auth/memory-adapter'
import PgAdapter from '@/app/extensions/next-auth/pg-adapter'

export default NextAuth({
  adapter: process.env.DATABASE_URL ? PgAdapter(process.env.DATABASE_URL) : MemoryAdapter(),
  providers: ([
    process.env.NEXTAUTH_GOOGLE ? GoogleProvider(JSON.parse(process.env.NEXTAUTH_GOOGLE)) : undefined,
    process.env.EMAIL_SERVER && process.env.EMAIL_FROM ? EmailProvider({
      server: process.env.EMAIL_SERVER,
      from: process.env.EMAIL_FROM,
    }) : undefined,
  ] as Array<Provider | undefined>).filter((v): v is Provider => v !== undefined),
})
