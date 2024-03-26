import React from 'react'
import { GPTAssistantCreate } from "@/app/api/client"
import { useAPIMutation } from "@/core/api/client"
import { useRouter } from "next/router"
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'
import * as Auth from 'next-auth/react'
import classNames from 'classnames'
import dynamic from 'next/dynamic'
import Head from 'next/head'
const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))

export default function ChatThread() {
  const router = useRouter()
  const auth = useSessionWithId()
  const {trigger, error, isMutating} = useAPIMutation(GPTAssistantCreate)
  return (
    <Layout>
      <Head><title>Chat</title></Head>
      <main className="flex-grow container mx-auto py-4 flex flex-col">
        {!auth.data?.user?.id ? (
          <div className="alert alert-warning shadow-lg block prose max-w-none">
            You are required to &nbsp; <button className="btn btn-sm" onClick={() => {Auth.signIn()}}>sign in</button> &nbsp; to use the chat.
          </div>
        ) : (
          <div className="prose">
            <h1>Chat</h1>
            <p>This is an experimental chat interface that translates your instructions to PWB workflows. Click start to begin!</p>
            <div className="flex flex-row items-center">
              <div className="btn btn-lg btn-primary" onClick={async (evt) => {
                const thread_id = await trigger()
                router.push(`/chat/${thread_id}`)
              }}>Start</div>
            </div>
            {error ? <div className="alert alert-error">{error.toString()}</div> : null}
            <progress className={classNames('progress progress-primary', { 'hidden': !isMutating })} />
          </div>
        )}
      </main>
    </Layout>
  )
}