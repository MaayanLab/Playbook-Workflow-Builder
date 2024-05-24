import React from 'react'
import { GPTAssistantDelete, GPTAssistantMessage, GPTAssistantMessagesList } from "@/app/api/client"
import { useAPIMutation, useAPIQuery } from "@/core/api/client"
import dynamic from "next/dynamic"
import Head from "next/head"
import * as Auth from 'next-auth/react'
import { useRouter } from "next/router"
import classNames from 'classnames'
import krg from '@/app/krg'
import * as dict from '@/utils/dict'
import type { Session } from 'next-auth'
import { AssembleState } from '@/app/api/v1/chat/utils'
import SafeRender from '@/utils/saferender'

import type { ReactMarkdownOptions } from 'react-markdown/lib/react-markdown'
const ReactMarkdown = dynamic(() => import('react-markdown/lib/react-markdown').then(({ ReactMarkdown }) => ReactMarkdown as ((props: ReactMarkdownOptions) => React.ReactNode)), { ssr: false })
const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const UserAvatar = dynamic(() => import('@/app/fragments/playbook/avatar'))

type Component = {
  id: number,
  inputs: Record<string, { id: number }>,
  name: string,
  value?: any,
}

function Component({ state, setState, component }: {
  state: Record<number, any>, setState: React.Dispatch<React.SetStateAction<Record<number, any>>>,
  component: Component,
}) {
  const metanode = React.useMemo(() => krg.getProcessNode(component.name), [component])
  const inputs = React.useMemo(() => {
    if (dict.keys(metanode.inputs).some(arg => !(`${component.inputs[arg].id}` in state))) {
      return undefined
    }
    const inputs: Record<string, unknown> = {}
    dict.items(metanode.inputs).forEach(({ key, value }) => {
      inputs[key] = value.codec.decode(state[`${component.inputs[key].id}`])
    })
    return inputs
  }, [metanode, state])
  const output = React.useMemo(() => {
    let output
    try{
      output = typeof state[`${component.id}`] !== 'undefined' ? metanode.output.codec.decode(state[`${component.id}`]) : component.value
    } catch (e) {}
    return output
  }, [state, component])
  React.useEffect(() => {
    if (`${component.id}` in state) return
    if (component.value !== undefined) {
      setState(state => ({...state, [`${component.id}`]: metanode.output.codec.encode(component.value)}))
      return
    }
    if (!('resolve' in metanode)) return
    if (!inputs) return
    (async () => {
      let result: any
      try {
        const req = await fetch(`/api/resolver`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            spec: metanode.spec,
            data: state[`${component.id}`],
            inputs,
          }),
        })
        if (!req.ok) throw new Error('Error')
        result = await req.json()
      } catch (e) { console.warn(e); result = '' }
      setState(state => ({ ...state, [`${component.id}`]: result }))
    })()
  }, [state, component])
  if (!metanode) return <span>Invalid metanode {component.name}</span>
  if ('prompt' in metanode) {
    return <div className="flex flex-col mx-16">
      <div className="card card-bordered rounded-3xl p-4">
        <div className="card-title flex-col place-items-start">
          <h3 className="prose text-xl">{metanode.meta.label}</h3>
          <h5 className="prose text-md">{metanode.story({ inputs, output })}</h5>
        </div>
        <div className="card-body">
          {inputs ?
            <SafeRender
              component={metanode.prompt}
              props={{
                inputs,
                output,
                submit: (output) => {
                  setState(state => ({ [`${component.id}`]: metanode.output.codec.encode(output) }))
                },
              }}
            />
            : <>Waiting for inputs...</>}
        </div>
      </div>
      {/* <div className="chat-footer text-lg cursor-pointer">
        <div className="bp5-icon bp5-icon-link" /> <div className="bp5-icon bp5-icon-thumbs-up" /> <div className="bp5-icon bp5-icon-thumbs-down" />
      </div> */}
    </div>
  }
  return <div className="flex flex-col mx-16">
    <div className="card card-bordered rounded-3xl p-4">
      <div className="card-title flex-col place-items-start">
        <h3 className="prose text-xl">{metanode.meta.label}</h3>
        <h5 className="prose text-md">{metanode.story({ inputs, output })}</h5>
      </div>
      <div className="card-body">
        {output ? <SafeRender component={metanode.output.view} props={metanode.output.codec.decode(output)} /> : <>Waiting...</>}
      </div>
    </div>
    {/* <div className="chat-footer text-lg cursor-pointer">
      <div className="bp5-icon bp5-icon-link" /> <div className="bp5-icon bp5-icon-thumbs-up" /> <div className="bp5-icon bp5-icon-thumbs-down" />
    </div> */}
  </div>
}

function Message({ session, role, children }: React.PropsWithChildren<{ session: Session | null, role: string }>) {
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
        {/* {role !== 'user' ?
          <div className={classNames('chat-footer text-lg cursor-pointer')}>
            <div className="bp5-icon bp5-icon-link" /> <div className="bp5-icon bp5-icon-thumbs-up" /> <div className="bp5-icon bp5-icon-thumbs-down" />
          </div>
          : null} */}
      </div>
    </>
  )
}

export default function ChatThread() {
  const { data: session } = Auth.useSession()
  const router = useRouter()
  const [state, setState] = React.useState({} as Record<number, string>)
  const [message, setMessage] = React.useState('')
  const { data: messages, mutate } = useAPIQuery(GPTAssistantMessagesList, { thread_id: router.query.thread_id as string })
  const { trigger, isMutating } = useAPIMutation(GPTAssistantMessage, { thread_id: router.query.thread_id as string })
  const { trigger: triggerDelete } = useAPIMutation(GPTAssistantDelete, { thread_id: router.query.thread_id as string })
  const playbookState = React.useMemo(() => messages ? AssembleState(messages, { with_value: true }) : undefined, [messages])
  const submit = React.useCallback(async (body: { message: string } | { step: { id: number, value?: string } }) => {
    const newMessages = await trigger({ body })
    await mutate(messages => [...messages ?? [], ...newMessages ?? []])
  }, [trigger])
  return (
    <Layout>
      <Head><title>Chat</title></Head>
      <main className="flex-grow container mx-auto p-4 flex flex-col gap-6">
        <div className="flex-grow prose max-w-none flex flex-row justify-between">
          <h1>Chat</h1>
          {/* <button
            type="button"
            className="btn btn-sm btn-error"
            onClick={() => {triggerDelete().finally(() => {router.push('/chat')})}}
          >Close</button> */}
        </div>
        <Message role="welcome" session={session}>
          I'm an AI-powered chat assistant interface designed to help you access the functionality of the playbook workflow builder.
          Please start by asking your question of interest, and I'll try my best to help you answer it through the construction of a playbook workflow.
        </Message>
        <div className="flex flex-row flex-wrap justify-center gap-2 place-self-center">
          {[
            'Show me the expression of ACE2 in healthy human tissues from GTEx',
            'Find drugs from the LINCS L1000 Chemical Perturbations the up regulate STAT3',
          ].map((suggestion, i) => {
            return (
              <button
                key={i}
                className="btn btn-ghost border border-primary btn-rounded rounded-lg btn-sm"
                onClick={evt => {submit({ message: suggestion })}}
              >{suggestion}</button>
            )
          })}
        </div>
        {messages?.map((message, i) => {
          const component = 'step' in message && playbookState ? playbookState.all_nodes[message.step.id] : undefined
          return (
            <React.Fragment key={i}>
              {'message' in message ? <Message role={message.role} session={session}>{message.message}</Message> : null}
              {component ? <Component component={component} state={state} setState={setState} /> : null}
              {message.role === 'assistant' && message.suggestions.length > 1 ?
                <div className="flex flex-row flex-wrap justify-center gap-2 place-self-center">
                  {message.suggestions.map((suggestion: any) => {
                    const suggestionProcess = playbookState?.all_nodes[suggestion.id]
                    const suggestionNode = suggestionProcess ? krg.getProcessNode(suggestionProcess.name) : undefined
                    if (!suggestionNode) return null
                    return (
                      <div key={suggestion.id} className="tooltip" data-tip={suggestionNode.meta.description}>
                        <button
                          className="btn btn-ghost border border-primary btn-rounded rounded-lg btn-sm"
                          onClick={evt => {submit({ step: suggestion })}}
                        >{suggestionNode.meta.label}</button>
                      </div>
                    )
                  })}
                </div>
                : null}
            </React.Fragment>
          )
        })}
        {isMutating ?
          <Message role="assistant" session={session}>
            <progress className="progress w-full"></progress>
          </Message>
          : null}
        <Message role="user" session={session}>
          <form
            className="flex flex-row items-center"
            onSubmit={async (evt) => {
              evt.preventDefault()
              await submit({ message })
              setMessage(() => '')
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
