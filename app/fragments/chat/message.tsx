import React from "react"
import { GPTAssistantMessageFeedback } from "@/app/api/client"
import { useAPIMutation } from "@/core/api/client"
import classNames from "classnames"
import dynamic from "next/dynamic"
import type { ReactMarkdownOptions } from 'react-markdown/lib/react-markdown'
import type { Session } from 'next-auth'

const ReactMarkdown = dynamic(() => import('react-markdown/lib/react-markdown').then(({ ReactMarkdown }) => ReactMarkdown as ((props: ReactMarkdownOptions) => React.ReactNode)), { ssr: false })
const UserAvatar = dynamic(() => import('@/app/fragments/playbook/avatar'))

export default function Message({ thread_id, message_id, session, role, children }: React.PropsWithChildren<{ thread_id?: string, message_id?: string, session: Session | null, role: string }>) {
  const [feedback, setFeedback] = React.useState('')
  const { trigger } = useAPIMutation(GPTAssistantMessageFeedback, { thread_id, message_id })
  return (
    <>
      <div className={classNames('chat', { 'chat-end': role === 'user', 'chat-start': role !== 'user', 'hidden': role === 'system' })}>
        <div className="chat-image btn btn-ghost btn-circle avatar placeholder">
          <div className="bg-neutral-focus text-neutral-content rounded-full w-16">
            {role !== 'user' ? "PWB" : <UserAvatar session={session} />}
          </div>
        </div>
        <div className={classNames('chat-bubble rounded-xl w-full prose', { 'chat-bubble-primary': role === 'user', 'chat-bubble-secondary': role === 'assistant' || role === 'welcome', 'chat-bubble-error': role === 'error' })}>
          {typeof children === 'string' ? <ReactMarkdown>{children}</ReactMarkdown> : children}
        </div>
        {role === 'assistant' ?
          <div className={classNames('chat-footer text-lg cursor-pointer')}>
            {/* <div
              className={classNames('bp5-icon bp5-icon-link', {'text-green-500': feedback === '+1'})}
            /> */}
            &nbsp;
            <div
              className={classNames('bp5-icon bp5-icon-thumbs-up', {'text-green-500': feedback === '+1'})}
              onClick={() => {
                trigger({ body: '+1' }).then(() => setFeedback('+1'))
              }}
            />
            &nbsp;
            <div
              className={classNames('bp5-icon bp5-icon-thumbs-down', {'text-red-500': feedback === '-1'})}
              onClick={() => {
                trigger({ body: '-1' }).then(() => setFeedback('-1'))
              }}
            />
          </div>
          : null}
      </div>
    </>
  )
}
