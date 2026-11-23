import React from "react"
import { useMutation } from '@tanstack/react-query'
import { useTRPC } from '@/app/trpc/client'
import classNames from "classnames"
import dynamic from "next/dynamic"
import type { Options as ReactMarkdownOptions } from 'react-markdown'
import type { Session } from '@/app/extensions/auth-client'
import { ExLink } from "@/app/fragments/ex-router"

const ReactMarkdown = dynamic<ReactMarkdownOptions>(() => import('react-markdown'), { ssr: false })
/** Tool call arguments/outputs are stored as JSON strings (older UI updates stored objects) */
function parseField(value: unknown): any {
  if (typeof value === 'string') return JSON.parse(value)
  return value ?? {}
}

const UserDisplay = dynamic(() => import('@/app/fragments/playbook/avatar').then(({ UserDisplay }) => UserDisplay))

export default function Message({ thread_id, message_id, session, role, children }: React.PropsWithChildren<{ thread_id?: string, message_id?: string, session: Session | null, role: string, embedded?: boolean }>) {
  const [feedback, setFeedback] = React.useState('')
  const trpc = useTRPC()
  const { mutateAsync: sendFeedback } = useMutation(trpc.chat.feedback.mutationOptions())
  if (role === 'developer' && typeof children === 'string') {
    try {
      const data = JSON.parse(children)
      if (data.function_call.name === 'options') {
        const args = parseField(data.function_call.arguments)
        if (args.workflow_id) {
          return <div className="p-1 my-1 italic">Checked possibilities from <ExLink href={`/graph/${args.workflow_id}/expand`}>here</ExLink></div>
        } else {
          return <div className="p-1 my-1 italic">Checking input possibilities</div>
        }
      } else if (data.function_call.name === 'search_published') {
        const args = parseField(data.function_call.arguments)
        if (args.search) {
          return <div className="p-1 my-1 italic">Searched for published workflows about "{args?.search}"</div>
        } else {
          return <div className="p-1 my-1 italic">Searching for published workflows</div>
        }
      } else if (data.function_call.name === 'expand') {
        const output = parseField(data.function_call_output?.output)
        if (output.error) {
          return <div className="p-1 my-1 italic">Agent experienced an error trying to extend workflow</div>
        } else if (output.dry_run) {
          return <div className="p-1 my-1 italic">Checked a step to add to the workflow</div>
        } else if (output.result?.label || output.result?.type) {
          return <div className="p-1 my-1 italic">Added {output.result.label ?? 'step'} to <ExLink href={`/report/${output.result.workflow_id}`}>the workflow</ExLink>{output.result.type ? ` to get ${output.result.type}` : ''}</div>
        } else {
          return <div className="p-1 my-1 italic">Adding step to workflow</div>
        }
      } else if (data.function_call.name === 'plan') {
        const output = parseField(data.function_call_output?.output)
        if (output.error) {
          return <div className="p-1 my-1 italic">Agent experienced an error trying to extend workflow</div>
        } else if (output.dry_run) {
          return <div className="p-1 my-1 italic">Checked {output.dry_run.steps?.length ?? 'several'} steps to add to the workflow</div>
        } else {
          return <div className="p-1 my-1 italic">Added {output.result?.steps?.length ?? 'several'} steps to the workflow</div>
        }
      } else if (data.function_call.name === 'set_value') {
        const output = parseField(data.function_call_output?.output)
        if (output.error) {
          return <div className="p-1 my-1 italic">Agent experienced an error trying to change an input</div>
        } else {
          return <div className="p-1 my-1 italic">Changed the {output.result?.label ?? 'input'} of the workflow</div>
        }
      } else if (data.function_call.name === 'find_steps' || data.function_call.name === 'describe' || data.function_call.name === 'help') {
        return <div className="p-1 my-1 italic">Looked up available capabilities</div>
      } else if (data.function_call.name === 'extend') {
        // a step added in the UI, as recorded by older threads (these are now recorded as `expand`)
        const args = parseField(data.function_call.arguments)
        const output = parseField(data.function_call_output?.output)
        const workflow_id = output.result?.workflow_id
        const what = args.label ?? args.step?.type ?? 'a step'
        return <div className="p-1 my-1 italic">Added {what} to {workflow_id ? <ExLink href={`/report/${workflow_id}`}>the workflow</ExLink> : 'the workflow'}</div>
      } else if (data.function_call.name === 'view') {
        const args = parseField(data.function_call.arguments)
        if (args.workflow_id) {
          return <div className="p-1 my-1 italic">Viewing <ExLink href={`/report/${args.workflow_id}`}>the workflow</ExLink></div>
        } else {
          return <div className="p-1 my-1 italic">Attempting to view workflow</div>
        }
      } else if (typeof data.function_call?.name === 'string') {
        // a tool without a description here: never show its raw JSON
        return <div className="p-1 my-1 italic">Agent used {data.function_call.name}</div>
      }
    } catch (error: any) {
      return <div className="p-1 my-1 italic">Agent produced a malformed request</div>
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
              className={classNames('bp6-icon bp6-icon-link', {'text-green-500': feedback === '+1'})}
            /> */}
            &nbsp;
            <div
              className={classNames('bp6-icon bp6-icon-thumbs-up', {'text-green-500': feedback === '+1'})}
              onClick={() => {
                if (thread_id && message_id) sendFeedback({ thread_id, message_id, feedback: '+1' }).then(() => setFeedback('+1'))
              }}
            />
            &nbsp;
            <div
              className={classNames('bp6-icon bp6-icon-thumbs-down', {'text-red-500': feedback === '-1'})}
              onClick={() => {
                if (thread_id && message_id) sendFeedback({ thread_id, message_id, feedback: '-1' }).then(() => setFeedback('-1'))
              }}
            />
          </div>
          : null}
      </div>
    </>
  )
}
