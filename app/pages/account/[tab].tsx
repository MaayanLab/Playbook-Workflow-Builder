import { useSessionWithId } from '@/app/extensions/next-auth/hooks'
import dynamic from 'next/dynamic'
import Head from 'next/head'
import React from 'react'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const AccountLayout = dynamic(() => import('@/app/fragments/account/layout'))

export default function Account() {
  const { data: session } = useSessionWithId({ required: true })
  return (
    <Layout>
      <Head>
        <title>Playbook: Account</title>
      </Head>

      <main className="flex-grow container mx-auto py-4 flex flex-row">
        {(session && session.user) ?
          <AccountLayout session={session} />
          : null}
      </main>
    </Layout>
  )
}
