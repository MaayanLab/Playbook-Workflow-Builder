import React from "react"
import { GPTAssistantMessageFeedback } from "@/app/api/client"
import { useAPIMutation } from "@/core/api/client"
import classNames from "classnames"
import dynamic from "next/dynamic"
import type { ReactMarkdownOptions } from 'react-markdown/lib/react-markdown'
import type { Session } from 'next-auth'

const ReactMarkdown = dynamic(() => import('react-markdown/lib/react-markdown').then(({ ReactMarkdown }) => ReactMarkdown as ((props: ReactMarkdownOptions) => React.ReactNode)), { ssr: false })
const UserDisplay = dynamic(() => import('@/app/fragments/playbook/avatar').then(({ UserDisplay }) => UserDisplay))

export default function Message({ thread_id, message_id, session, role, children, embedded }: React.PropsWithChildren<{ thread_id?: string, message_id?: string, session: Session | null, role: string, embedded?: boolean }>) {
  const [feedback, setFeedback] = React.useState('')
  const { trigger } = useAPIMutation(GPTAssistantMessageFeedback, { thread_id, message_id })
  if (embedded) {
    if (role === 'developer' && typeof children === 'string') {
      try {
        const data = JSON.parse(children)
        if (data.function_call.name === 'extend') {
          const output = JSON.parse(data.function_call_output.output)
          return <div>Added step to workflow to get {output.result.type}</div>
        } else {
          return null
        }
      } catch (error: any) {
        console.error({ children, error })
        return null
      }
    } else if (role === 'reasoning' && typeof children === 'string') {
      return <div><i>{children}</i></div>
    }
  }
  return (
    <>
      <div className={classNames('chat', { 'chat-end': role === 'user', 'chat-start': role !== 'user' })}>
        <div className="chat-header">
          {role !== 'user' ? "PWB" : <UserDisplay session={session} />}
        </div>
        <div className={classNames('chat-bubble w-full prose', { 'chat-bubble-primary': role === 'user', 'chat-bubble-secondary': role === 'assistant' || role === 'welcome', 'chat-bubble-error': role === 'error' })}>
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
