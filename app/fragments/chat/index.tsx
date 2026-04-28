import React from 'react'
import dynamic from "next/dynamic"
import Head from "next/head"
import fetcher from '@/utils/next-rest-fetcher'
import { MetapathProvider } from '@/app/fragments/metapath'
import { SWRConfig } from 'swr'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const Chat = dynamic(() => import('@/app/fragments/chat/chat'))

export default function ChatThread() {
  return (
    <Layout>
      <Head><title>Text to Workflow</title></Head>
      <SWRConfig value={{ fetcher }}>
        <MetapathProvider>
          <main className="flex-grow container mx-auto p-4 flex flex-col gap-6">
            <Chat mode="report" />
          </main>
        </MetapathProvider>
      </SWRConfig>
    </Layout>
  )
}
