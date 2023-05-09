import React from 'react'
import dynamic from 'next/dynamic'
import Head from 'next/head'
import * as Auth from 'next-auth/react'
import classNames from 'classnames'
import { ReactMarkdown } from 'react-markdown/lib/react-markdown'
import krg from '@/app/krg'
import * as dict from '@/utils/dict'
import type { Session } from 'next-auth'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const UserAvatar = dynamic(() => import('@/app/fragments/playbook/avatar'))

function Workflow({ workflow }: { workflow: { spec: string, data?: string }[] }) {
  const [state, setState] = React.useState({} as Record<number, string>)
  React.useEffect(() => {
    setState(() => dict.init(
      workflow
        .filter((workflow): workflow is { spec: string, data: string } => workflow.data !== undefined)
        .map(({ data }, i) => ({ key: i, value: data }))))
  }, [workflow])
  React.useEffect(() => {
    if (Object.keys(state).length > 0 && Object.keys(state).length < workflow.length) {
      console.log({ state, workflow })
      const n = Math.max(...Object.keys(state).map(k => +k))
      const proc = krg.getProcessNode(workflow[n+1].spec)
      if ('resolve' in proc) {
        const formData = new FormData()
        dict.items(proc.inputs).forEach(({ key, value }) => {
          formData.append(key, value.codec.encode(state[n]))
        })
        const req = fetch(`/api/resolver/${proc.spec}`, {
          method: 'POST',
          body: formData,
        }).then(req => req.json())
        .then(res => {
          let result: any
          try { result = proc.output.codec.decode(res) }
          catch (e) { result = '' }
          setState(state => ({ ...state, [n+1]: result }))
        })
      }
    }
  }, [state, workflow.length])
  return (
    <div className='flex flex-col gap-2'>
      {workflow.map(({ spec, data }, i) => {
        const metanode = krg.getProcessNode(spec)
        if (!metanode) return <span>Invalid metanode {spec}</span>
        if ('prompt' in metanode) {
          const Prompt = metanode.prompt
          return <div key={i}>
            <h3 className="m-0">{metanode.meta.label}</h3>
            <Prompt inputs={{}} output={state[i]} submit={(output) => {
              setState((data) => {
                const newData = {} as Record<number, string>
                for (const k in data) {
                  if (+k >= i) continue
                  newData[k] = data[k]
                }
                newData[i] = output as string
                return newData
              })
            }} />
          </div>
        } else {
          let viewOutput
          try {
            viewOutput = metanode.output.view(state[i])
          } catch (e) {
            viewOutput = <div>Error: {e.toString()}</div>
          }
          return <div key={i}>
            <h3 className="m-0">{metanode.meta.label}</h3>
            {viewOutput}
          </div>
        }
      })}
    </div>
  )
}

function Message({ session, role, children }: React.PropsWithChildren<{ session: Session | null, role: 'user' | 'assistant' | 'system' | 'error' }>) {
  return (
    <div className={classNames('chat', { 'chat-end': role === 'user', 'chat-start': role !== 'user', 'hidden': role === 'system' })}>
      <div className="chat-image btn btn-ghost btn-circle avatar placeholder">
        <div className="bg-neutral-focus text-neutral-content rounded-full w-16">
          {role !== 'user' ? "PWB" : <UserAvatar session={session} />}
        </div>
      </div>
      <div className={classNames('chat-bubble rounded-xl w-full prose', { 'chat-bubble-primary': role === 'user', 'chat-bubble-secondary': role === 'assistant', 'chat-bubble-error': role === 'error' })}>
        {typeof children === 'string' ?
          <ReactMarkdown
            components={React.useMemo(() => ({
              pre({ className, children, ...props }) {
                return <div {...props} className={classNames(className, 'bg-transparent')}>{children}</div>
              },
              code({ node, inline, className, children, ...props }) {
                const match = /language-(\w+)/.exec(className || '')
                if (!inline && match && match[1] === 'workflow') {
                  const workflow = JSON.parse(children.join('\n'))
                  return <Workflow workflow={workflow} />
                }
                return <code {...props} className={classNames(className, 'bg-secondary-content')}>{children}</code>
              }
            }), [])}
          >{children}</ReactMarkdown>
          : children}
      </div>
    </div>
  )
}

export default function Chat() {
  const { data: session } = Auth.useSession()
  const [message, setMessage] = React.useState('')
  const [chat, setChat] = React.useState({
    waitingForReply: false,
    messages: [] as { role: 'system' | 'user' | 'assistant' | 'error', content: string }[],
  })
  return (
    <Layout>
      <Head><title>Chat</title></Head>
      <main className="flex-grow container mx-auto p-4 flex flex-col gap-6 justify-end">
        <Message role="assistant" session={session}>
          I'm an AI-powered chat assistant interface designed to help you access the functionality of the playbook workflow builder. Please start by asking your question of interest, and I'll try my best to help you answer it through the construction of a playbook workflow.
        </Message>
        {chat.messages.map((message, i) =>
          <Message key={i} role={message.role} session={session}>{message.content}</Message>
        )}
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
                  messages: [...chat.messages, { role: 'user', content: message}].filter(({ role }) => role !== 'error'),
                }),
              })
              const res = req.ok ? await req.json() : { role: 'error', content: await req.text() }
              setChat((cc) => ({ waitingForReply: false, messages: [...cc.messages, res] }))
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
