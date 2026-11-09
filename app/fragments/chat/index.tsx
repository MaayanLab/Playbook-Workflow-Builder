import React from 'react'
import dynamic from "next/dynamic"
import DocumentTitle from '@/app/fragments/document-title'
import { MetapathProvider } from '@/app/fragments/metapath'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const Chat = dynamic(() => import('@/app/fragments/chat/chat'))

export default function ChatThread() {
  return (
    <>
      <MetapathProvider>
      <DocumentTitle title={'Text to Workflow'} />
        <Layout>
          <main className="flex-grow container mx-auto p-4 flex flex-col gap-6">
            <Chat mode="chat" />
          </main>
        </Layout>
      </MetapathProvider>
    </>
  )
}
