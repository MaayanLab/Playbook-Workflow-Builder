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

function Component({ state, setState, component }: {
  state: Record<number, string>, setState: React.Dispatch<React.SetStateAction<Record<number, string>>>,
  component: { id: number, inputs: Record<string, number>, type: string, description: string, data?: string },
}) {
  const metanode = krg.getProcessNode(component.type)
  React.useEffect(() => {
    if (component.id in state) return
    if (component.data !== undefined) {
      setState({[component.id]: metanode.output.codec.encode(component.data)})
      return
    }
    if (!('resolve' in metanode)) return
    if (dict.keys(metanode.inputs).some(arg => !(component.inputs[arg] in state))) return
    const formData = new FormData()
    dict.keys(metanode.inputs).forEach((arg) => {
      formData.append(arg, state[component.inputs[arg]])
    })
    fetch(`/api/resolver/${metanode.spec}`, {
      method: 'POST',
      body: formData,
    }).then(req => req.json())
    .then(res => {
      let result: any
      try { result = res }
      catch (e) { console.warn(e); result = '' }
      setState(state => ({ ...state, [component.id]: result }))
    })
  }, [state, component])
    if (!metanode) return <span>Invalid metanode {component.type}</span>
  if ('prompt' in metanode) {
    const Prompt = metanode.prompt
    const inputs: Record<string, unknown> = {}
    dict.items(metanode.inputs).forEach(({ key, value }) => {
      inputs[key] = value.codec.decode(state[component.inputs[key]])
    })
    const output = typeof state[component.id] !== undefined ? metanode.output.codec.decode(state[component.id]) : component.data
    return <div>
      <h3 className="m-0">{metanode.meta.label}</h3>
      <Prompt
        inputs={inputs}
        output={output}
        submit={(output) => {
          setState({ [component.id]: metanode.output.codec.encode(output) })
        }}
      />
    </div>
  }
  const output = state[component.id]
  let viewOutput
  try {
    viewOutput = output ? metanode.output.view(metanode.output.codec.decode(output)) : <div>Updating...</div>
  } catch (e) {
    viewOutput = <div>Error: {(e as any).toString()}</div>
  }
  return <div>
    <h3 className="m-0">{metanode.meta.label}</h3>
    {viewOutput}
  </div>
}

function MessageOutput({ children, state, setState }: React.PropsWithChildren<{
  state: Record<number, string>, setState: React.Dispatch<React.SetStateAction<Record<number, string>>>,
}>) {
  const { component } = React.useMemo(() => {
    if (typeof children === 'string') {
      const m = /```component\n(.+)\n```/g.exec(children)
      if (m) {
        return { component: JSON.parse(m[1]) }
      }
    }
    return {}
  }, [children])
  if (component) return <Component component={component} state={state} setState={setState} />
  else if (typeof children === 'string') return <ReactMarkdown>{children}</ReactMarkdown>
  else return <>{children}</>
}

function Message({
  session, role, children,
  state, setState,
}: React.PropsWithChildren<{
  session: Session | null, role: 'welcome' | 'user' | 'assistant' | 'system' | 'error'
  state: Record<number, string>, setState: React.Dispatch<React.SetStateAction<Record<number, string>>>,
}>) {
  return (
    <div className={classNames('chat', { 'chat-end': role === 'user', 'chat-start': role !== 'user', 'hidden': role === 'system' })}>
      <div className="chat-image btn btn-ghost btn-circle avatar placeholder">
        <div className="bg-neutral-focus text-neutral-content rounded-full w-16">
          {role !== 'user' ? "PWB" : <UserAvatar session={session} />}
        </div>
      </div>
      <div className={classNames('chat-bubble rounded-xl w-full prose', { 'chat-bubble-primary': role === 'user', 'chat-bubble-secondary': role === 'assistant' || role === 'welcome', 'chat-bubble-error': role === 'error' })}>
        <MessageOutput state={state} setState={setState}>{children}</MessageOutput>
      </div>
      {role === 'assistant' ? 
        <div className={classNames('chat-footer text-lg cursor-pointer')}>
          <div className="bp4-icon bp4-icon-link" /> <div className="bp4-icon bp4-icon-thumbs-up" /> <div className="bp4-icon bp4-icon-thumbs-down" />
        </div>
        : null}
    </div>
  )
}

export default function Chat() {
  const { data: session } = Auth.useSession()
  const [state, setState] = React.useState({} as Record<number, string>)
  const [message, setMessage] = React.useState('')
  const [chat, setChat] = React.useState({
    waitingForReply: false,
    messages: [] as { role: 'system' | 'user' | 'assistant' | 'error', content: string }[],
  })
  return (
    <Layout>
      <Head><title>Chat</title></Head>
      <main className="flex-grow container mx-auto p-4 flex flex-col gap-6 justify-end">
        <Message role="welcome" session={session} state={state} setState={setState}>
          I'm an AI-powered chat assistant interface designed to help you access the functionality of the playbook workflow builder. Please start by asking your question of interest, and I'll try my best to help you answer it through the construction of a playbook workflow.
        </Message>
        {chat.messages.map((message, i) =>
          <Message key={i} role={message.role} session={session} state={state} setState={setState}>{message.content}</Message>
        )}
        {chat.waitingForReply ?
          <Message role="assistant" session={session} state={state} setState={setState}>
            <progress className="progress w-full"></progress>
          </Message>
          : null}
        <Message role="user" session={session} state={state} setState={setState}>
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
