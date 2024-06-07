import React from 'react'
import { GPTAssistantCreate } from "@/app/api/client"
import { useAPIMutation } from "@/core/api/client"
import { useRouter } from "next/router"
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'
import * as Auth from 'next-auth/react'
import classNames from 'classnames'
import dynamic from 'next/dynamic'
import Head from 'next/head'
import usePublicUrl from '@/utils/next-public-url'
import Message from '@/app/fragments/chat/message'
const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))

export default function ChatThread() {
  const router = useRouter()
  const publicUrl = usePublicUrl()
  const auth = useSessionWithId()
  const {trigger, error, isMutating} = useAPIMutation(GPTAssistantCreate)
  React.useEffect(() => {
    if (!auth.data?.user?.id) return
    trigger().then((thread_id) => router.push(`/chat/${thread_id}`))
  }, [auth.data?.user?.id])
  return (
    <Layout>
      <Head><title>Chat</title></Head>
      <main className="flex-grow container mx-auto p-4 flex flex-col gap-6">
        <div className={classNames("flex-grow flex flex-col justify-center items-center")}>
          <img
            className="w-32"
            src={`${publicUrl}/PWB-logo.svg`}
          />
          <div className="prose"><h4>Playbook Workflow Builder Text to Workflow</h4></div>
          <div className="flex flex-col overflow-hidden self-stretch">
            <Message role="welcome" session={null}>
              How can I help you today?
            </Message>
            <progress className={classNames('progress progress-primary', { 'hidden': !isMutating })} />
            {error ? <div className="alert alert-error">{error.toString()}</div> : null}
            <Message role="user" session={null}>
              <div className="alert alert-warning shadow-lg block">
                You are required to &nbsp; <button className="btn btn-sm" onClick={() => {Auth.signIn()}}>sign in</button> &nbsp; to use the chat.
              </div>
            </Message>
          </div>
        </div>
      </main>
    </Layout>
  )
}