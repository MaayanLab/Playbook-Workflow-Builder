import React from 'react'
import dynamic from "next/dynamic"
import DocumentTitle from '@/app/fragments/document-title'

const Chat = dynamic(() => import('@/app/fragments/chat/chat'))

export default function ChatThread() {
  return (
    <>
      <DocumentTitle title={'Text to Workflow'} />
      <main className="flex-grow container mx-auto p-4 flex flex-col gap-6">
        <Chat mode="chat" />
      </main>
    </>
  )
}
