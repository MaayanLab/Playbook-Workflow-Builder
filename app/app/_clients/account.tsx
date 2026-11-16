'use client'

import { useSessionWithId } from '@/app/extensions/auth/hooks'
import dynamic from 'next/dynamic'
import React from 'react'

const AccountLayout = dynamic(() => import('@/app/fragments/account/layout'))

export default function Account() {
  const { data: session } = useSessionWithId({ required: true })
  return (
    <>

      <main className="flex-grow container mx-auto py-4 flex flex-row">
        {(session && session.user) ?
          <AccountLayout session={session} />
          : null}
      </main>
    </>
  )
}
