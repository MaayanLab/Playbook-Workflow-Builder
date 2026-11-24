import React from "react"
import { useMutation } from '@tanstack/react-query'
import { useTRPC } from '@/app/trpc/client'
import classNames from "classnames"
import dynamic from "next/dynamic"
import type { Options as ReactMarkdownOptions } from 'react-markdown'
import type { Session } from '@/app/extensions/auth-client'
import { ExLink } from "@/app/fragments/ex-router"

const ReactMarkdown = dynamic<ReactMarkdownOptions>(() => import('react-markdown'), { ssr: false })
/**
 * Tool call arguments/outputs are stored as JSON strings (older UI updates stored objects); while a
 *  call is in progress its arguments are empty (or partial) & it has no output yet.
 */
function parseField(value: unknown): any {
  if (typeof value === 'string') {
    try { return JSON.parse(value) } catch { return undefined }
  }
  return value ?? undefined
}

const Status = ({ children }: React.PropsWithChildren<{}>) => <div className="p-1 my-1 italic chat-message-in">{children}</div>
const Link = ({ workflow_id, children }: React.PropsWithChildren<{ workflow_id?: string }>) =>
  workflow_id ? <ExLink href={`/report/${workflow_id}`}>{children}</ExLink> : <>{children}</>

/** How a tool call is shown: in progress (no output yet), failed, or done */
function ToolCall({ name, args = {}, output, error }: { name: string, args?: any, output?: any, error?: unknown }) {
  const failed = Boolean(error || output?.error)
  const pending = !failed && output === undefined
  switch (name) {
    case 'options':
      if (failed) return <Status>Couldn't check the possibilities for the workflow</Status>
      return <Status>{pending ? 'Checking' : 'Checked'} {args.workflow_id ? <>possibilities from <ExLink href={`/graph/${args.workflow_id}/expand`}>here</ExLink></> : 'input possibilities'}</Status>
    case 'search_published':
      if (failed) return <Status>Couldn't search published workflows</Status>
      return <Status>{pending ? 'Searching' : 'Searched'} for published workflows{args.search ? <> about "{args.search}"</> : null}</Status>
    case 'expand':
      if (failed) return <Status>Agent experienced an error trying to extend the workflow</Status>
      if (pending) return <Status>Adding a step to the workflow</Status>
      if (output.dry_run) return <Status>Checked a step to add to the workflow</Status>
      return <Status>Added {output.result?.label ?? 'a step'} to <Link workflow_id={output.result?.workflow_id}>the workflow</Link>{output.result?.type ? ` to get ${output.result.type}` : ''}</Status>
    case 'plan':
      if (failed) return <Status>Agent experienced an error trying to extend the workflow</Status>
      if (pending) return <Status>Adding {args.steps?.length ?? 'several'} steps to the workflow</Status>
      if (output.dry_run) return <Status>Checked {output.dry_run.steps?.length ?? 'several'} steps to add to the workflow</Status>
      return <Status>Added {output.result?.steps?.length ?? 'several'} steps to <Link workflow_id={output.result?.workflow_id}>the workflow</Link></Status>
    case 'set_value':
      if (failed) return <Status>Agent experienced an error trying to change an input</Status>
      if (pending) return <Status>Changing an input of the workflow</Status>
      return <Status>Changed the {output.result?.label ?? 'input'} of <Link workflow_id={output.result?.workflow_id}>the workflow</Link></Status>
    case 'find_steps':
    case 'describe':
    case 'help':
      return <Status>{pending ? 'Looking up' : 'Looked up'} available capabilities</Status>
    case 'view':
      if (failed) return <Status>Couldn't view the workflow</Status>
      return <Status>{pending ? 'Viewing' : 'Viewed'} <Link workflow_id={args.workflow_id}>the workflow</Link></Status>
    case 'extend': {
      // a step added in the UI, as recorded by older threads (these are now recorded as `expand`)
      return <Status>Added {args.label ?? args.step?.type ?? 'a step'} to <Link workflow_id={output?.result?.workflow_id}>the workflow</Link></Status>
    }
    default:
      // a tool without a description here: never show its raw JSON
      return <Status>{pending ? 'Using' : 'Used'} {name}</Status>
  }
}

const UserDisplay = dynamic(() => import('@/app/fragments/playbook/avatar').then(({ UserDisplay }) => UserDisplay))

export default function Message({ thread_id, message_id, session, role, children }: React.PropsWithChildren<{ thread_id?: string, message_id?: string, session: Session | null, role: string, embedded?: boolean }>) {
  const [feedback, setFeedback] = React.useState('')
  const trpc = useTRPC()
  const { mutateAsync: sendFeedback } = useMutation(trpc.chat.feedback.mutationOptions())
  if (role === 'developer' && typeof children === 'string') {
    let data: any
    try { data = JSON.parse(children) } catch {}
    if (typeof data?.function_call?.name !== 'string') return <Status>Agent produced a malformed request</Status>
    return <ToolCall
      name={data.function_call.name}
      args={parseField(data.function_call.arguments)}
      output={parseField(data.function_call_output?.output)}
      error={data.function_call_output?.error}
    />
  }
  return (
    <>
      <div className={classNames('chat chat-message-in', { 'chat-end': role === 'user', 'chat-start': role !== 'user' })}>
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
