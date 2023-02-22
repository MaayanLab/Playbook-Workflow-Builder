import { useSessionWithId } from '@/app/extensions/next-auth/hooks'
import dynamic from 'next/dynamic'
import Head from 'next/head'
import React from 'react'

const Header = dynamic(() => import('@/app/fragments/playbook/header'))
const Footer = dynamic(() => import('@/app/fragments/playbook/footer'))
const AccountLayout = dynamic(() => import('@/app/fragments/account/layout'))

export default function Account() {
  const { data: session } = useSessionWithId({ required: true })
  return (
    <div className="flex flex-col min-w-screen min-h-screen">
      <Head>
        <title>Playbook: Account</title>
      </Head>

      <Header />
      
      <main className="flex-grow container mx-auto py-4 flex flex-row">
        {(session && session.user) ?
          <AccountLayout session={session} />
          : null}
      </main>

      <Footer />
    </div>
  )
}
