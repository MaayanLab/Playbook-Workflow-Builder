import React from 'react'
import { GPTAssistantThreadDelete, GPTAssistantThreadMessage, GPTAssistantThreadMessagesList } from "@/app/api/client"
import { useAPIMutation, useAPIQuery } from "@/core/api/client"
import dynamic from "next/dynamic"
import Head from "next/head"
import * as Auth from 'next-auth/react'
import { useRouter } from "next/router"

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const Message = dynamic(() => import('@/app/fragments/chat/message'))

export default function ChatThread() {
  const { data: session } = Auth.useSession()
  const router = useRouter()
  const [state, setState] = React.useState({} as Record<number, string>)
  const [message, setMessage] = React.useState('')
  const { data: messages, mutate } = useAPIQuery(GPTAssistantThreadMessagesList, { thread_id: router.query.thread_id as string })
  const { trigger, isMutating } = useAPIMutation(GPTAssistantThreadMessage, { thread_id: router.query.thread_id as string })
  const { trigger: triggerDelete } = useAPIMutation(GPTAssistantThreadDelete, { thread_id: router.query.thread_id as string })
  const submit = React.useCallback(async ({ content }: { content: string }) => {
    const newMessages = await trigger({ body: { content } })
    mutate(messages => [...messages ?? [], ...newMessages?.data ?? []])
  }, [trigger])
  return (
    <Layout>
      <Head><title>Chat</title></Head>
      <main className="flex-grow container mx-auto p-4 flex flex-col gap-6">
        <div className="flex-grow text-right">
          <button
            type="button"
            className="btn btn-sm btn-error"
            onClick={() => {triggerDelete().finally(() => {router.push('/chat')})}}
          >Close</button>
        </div>
        <Message role="welcome" submit={submit} session={session} state={state} setState={setState}>
          I'm an AI-powered chat assistant interface designed to help you access the functionality of the playbook workflow builder.
          Please start by asking your question of interest, and I'll try my best to help you answer it through the construction of a playbook workflow.
        </Message>
        {messages?.map((message, i) => <Message key={i} role={message.role} submit={submit} session={session} state={state} setState={setState}>{message.content.map(content => content.type === 'text' ? <>{content.text.value}</> : null)}</Message>)}
        {isMutating ?
          <Message role="assistant" submit={submit} session={session} state={state} setState={setState}>
            <progress className="progress w-full"></progress>
          </Message>
          : null}
        <Message role="user" submit={submit} session={session} state={state} setState={setState}>
          <form
            className="flex flex-row items-center"
            onSubmit={async (evt) => {
              evt.preventDefault()
              submit({ content: message })
            }}
          >
            <input
              type="text"
              className="input w-full bg-transparent rounded-full"
              placeholder="Type your questions here"
              value={message}
              onChange={evt => setMessage(() => evt.target.value)}
            />
            <button type="submit" className="btn btn-sm" disabled={!message}>Send</button>
          </form>
        </Message>
      </main>
    </Layout>
  )
}
