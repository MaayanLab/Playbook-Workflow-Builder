import React from 'react'
import dynamic from 'next/dynamic'
import Head from 'next/head'
import * as Auth from 'next-auth/react'
import classNames from 'classnames'
import type { Session } from 'next-auth'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const UserAvatar = dynamic(() => import('@/app/fragments/playbook/avatar'))

function Message({ session, role, children }: React.PropsWithChildren<{ session: Session | null, role: 'user' | 'assistant' | 'error' }>) {
  return (
    <div className={classNames('chat', { 'chat-end': role === 'user', 'chat-start': role !== 'user' })}>
      <div className="chat-image btn btn-ghost btn-circle avatar placeholder">
        <div className="bg-neutral-focus text-neutral-content rounded-full w-16">
          {role !== 'user' ? "PWB" : <UserAvatar session={session} />}
        </div>
      </div>
      <div className={classNames('chat-bubble rounded-xl w-full prose', { 'chat-bubble-primary': role === 'user', 'chat-bubble-secondary': role === 'assistant', 'chat-bubble-error': role === 'error' })}>
        {children}
      </div>
    </div>
  )
}

export default function Chat() {
  const { data: session } = Auth.useSession()
  const [message, setMessage] = React.useState('')
  const [chat, setChat] = React.useState({
    waitingForReply: false,
    messages: [] as { role: 'user' | 'assistant' | 'error', content: string }[],
  })
  return (
    <Layout>
      <Head><title>Chat</title></Head>
      <main className="flex-grow container mx-auto p-4 flex flex-col gap-6 justify-end">
        <Message role="assistant" session={session}>
          I'm an AI-powered chat assistant interface designed to help you access the functionality of the playbook workflow builder. Please start by asking your question of interest, and I'll try my best to help you answer it through the construction of a playbook workflow.
        </Message>
        {chat.messages.map((message, i) => <Message key={i} role={message.role} session={session}>{message.content}</Message>)}
        {chat.waitingForReply ?
          <Message role="assistant" session={session}>
            <progress className="progress w-full"></progress>
          </Message>
          : null}
        <Message role="user" session={session}>
          <form
            className="flex flex-row items-center"
            onSubmit={async (evt) => {
              evt.preventDefault()
              if (chat.waitingForReply) return
              setChat((cc) => ({ waitingForReply: true, messages: [...cc.messages, { role: 'user', content: message }] }))
              setMessage(() => '')
              const req = await fetch('/api/gpt/chat', {
                method: 'POST',
                body: JSON.stringify({
                  messages: chat.messages.filter(({ role }) => role !== 'error'),
                }),
              })
              const msg = req.ok ? await req.json() : { role: 'error', content: await req.text() }
              setChat((cc) => ({ waitingForReply: false, messages: [...cc.messages, msg] }))
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
