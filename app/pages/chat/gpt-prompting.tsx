import React from 'react'
import dynamic from 'next/dynamic'
import Head from 'next/head'
import * as Auth from 'next-auth/react'
import { useAPIMutation } from '@/core/api/client'
import { ChatGPTPrompting } from '@/app/api/client'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const Message = dynamic(() => import('@/app/fragments/chat/message'))

export default function Chat() {
  const { data: session } = Auth.useSession()
  const [state, setState] = React.useState({} as Record<number, string>)
  const { trigger } = useAPIMutation(ChatGPTPrompting)
  const [message, setMessage] = React.useState('')
  const [chat, setChat] = React.useState({
    waitingForReply: false,
    messages: [] as { role: string, content: string }[],
  })
  const submit = React.useCallback(async (message: { role: string, content: string }) => {
    if (chat.waitingForReply) return
    setChat((cc) => ({ waitingForReply: true, messages: [...cc.messages, message] }))
    setMessage(() => '')
    const results = await trigger({
      body: {
        messages: [...chat.messages, message].filter(({ role }) => role !== 'error'),
      }
    })
    setChat((cc) => {
      if (!results) return { waitingForReply: false, messages: cc.messages }
      else return ({ waitingForReply: false, messages: [...cc.messages, ...results] })
    })
  }, [chat])
  return (
    <Layout>
      <Head><title>Chat</title></Head>
      <main className="flex-grow container mx-auto p-4 flex flex-col gap-6 justify-end">
        <Message role="welcome" submit={submit} session={session} state={state} setState={setState}>
          I'm an AI-powered chat assistant interface designed to help you access the functionality of the playbook workflow builder.
          Please start by asking your question of interest, and I'll try my best to help you answer it through the construction of a playbook workflow.
        </Message>
        {chat.messages.map((message, i) =>
          <Message key={i} role={message.role} submit={submit} session={session} state={state} setState={setState}>{message.content}</Message>
        )}
        {chat.waitingForReply ?
          <Message role="assistant" submit={submit} session={session} state={state} setState={setState}>
            <progress className="progress w-full"></progress>
          </Message>
          : null}
        <Message role="user" submit={submit} session={session} state={state} setState={setState}>
          <form
            className="flex flex-row items-center"
            onSubmit={async (evt) => {
              evt.preventDefault()
              submit({ role: 'user', content: message })
            }}
          >
            <input
              type="text"
              className="input w-full bg-transparent rounded-full"
              placeholder="Type your questions here"
              value={message}
              onChange={evt => setMessage(() => evt.target.value)}
            />
            <button type="submit" className="btn btn-sm" disabled={!message && !chat.waitingForReply}>Send</button>
          </form>
        </Message>
      </main>
    </Layout>
  )
}
