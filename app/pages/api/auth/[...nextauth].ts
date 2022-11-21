import NextAuth from 'next-auth'
import EmailProvider from 'next-auth/providers/email'
import kvdb from '@/app/kvdb'
import LevelAdapter from '@/app/extensions/next-auth/level-adapter'

export default NextAuth({
  adapter: LevelAdapter(kvdb),
  providers: [
    EmailProvider({
      server: process.env.EMAIL_SERVER,
      from: process.env.EMAIL_FROM,
    })
  ],
})