import React from 'react'
import { GPTAssistantThread2Delete, GPTAssistantThread2Message, GPTAssistantThread2MessagesList } from "@/app/api/client"
import { useAPIMutation, useAPIQuery } from "@/core/api/client"
import dynamic from "next/dynamic"
import Head from "next/head"
import * as Auth from 'next-auth/react'
import { useRouter } from "next/router"
import classNames from 'classnames'
import { ReactMarkdown } from 'react-markdown/lib/react-markdown'
import krg from '@/app/krg'
import * as dict from '@/utils/dict'
import type { Session } from 'next-auth'
import { z } from 'zod'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const UserAvatar = dynamic(() => import('@/app/fragments/playbook/avatar'))

type Component = {
  id: number,
  inputs: Record<string, { id: number }>,
  name: string,
  data?: any,
  likelihood?: number,
}

function Component({ state, setState, component }: {
  state: Record<number, any>, setState: React.Dispatch<React.SetStateAction<Record<number, any>>>,
  component: Component,
}) {
  const metanode = krg.getProcessNode(component.name)
  React.useEffect(() => {
    if (component.id in state) return
    if (component.data !== undefined) {
      setState({[component.id]: metanode.output.codec.encode(component.data)})
      return
    }
    if (!('resolve' in metanode)) return
    if (dict.keys(metanode.inputs).some(arg => !(component.inputs[arg].id in state))) return
    async () => {
      let result: any
      try {
        const req = await fetch(`/api/resolver`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            spec: metanode.spec,
            data: state[component.id],
            inputs: dict.init(dict.keys(metanode.inputs).map((arg) => ({ key: arg, value: state[component.inputs[arg].id] }))),
          }),
        })
        if (!req.ok) throw new Error('Error')
        result = await req.json()
      } catch (e) { console.warn(e); result = '' }
      setState(state => ({ ...state, [component.id]: result }))
    }
  }, [state, component])
    if (!metanode) return <span>Invalid metanode {component.name}</span>
  if ('prompt' in metanode) {
    const Prompt = metanode.prompt
    const inputs: Record<string, unknown> = {}
    dict.items(metanode.inputs).forEach(({ key, value }) => {
      inputs[key] = value.codec.decode(state[component.inputs[key].id])
    })
    let output
    try{
      output = typeof state[component.id] !== 'undefined' ? metanode.output.codec.decode(state[component.id]) : component.data
    } catch (e) {}
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
        {role === 'assistant' ? 
          <div className={classNames('chat-footer text-lg cursor-pointer')}>
            <div className="bp4-icon bp4-icon-link" /> <div className="bp4-icon bp4-icon-thumbs-up" /> <div className="bp4-icon bp4-icon-thumbs-down" />
          </div>
          : null}
      </div>
    </>
  )
}

function TryJSONParse(obj: string) {
  try {
    return JSON.parse(obj)
  } catch (e) {
    return obj
  }
}

function AssembleState(messages: (string | undefined)[]) {
  let max_id = 0
  const all_nodes: Record<number, { id: number, name: string, inputs: Record<string, { id: number }> }> = {}
  const all_values: Record<number, string> = {}
  const workflow: { id: number, name: string, inputs: Record<string, { id: number }> }[] = []
  messages
    .filter((msg): msg is string => !!msg)
    .map(TryJSONParse)
    .map(msg => z.union([
      z.string().transform(el => ({ message: el })),
      z.object({ message: z.string(), choices: z.array(z.object({ id: z.number(), name: z.string(), inputs: z.record(z.string(), z.object({ id: z.number() })) })) }),
      z.object({ step: z.object({ id: z.number() }), choices: z.array(z.object({ id: z.number(), name: z.string(), inputs: z.record(z.string(), z.object({ id: z.number() })) })) }),
    ]).parse(msg))
    .forEach(item => {
      if ('step' in item) {
        workflow.push(all_nodes[+item.step.id])
      }
      if ('choices' in item) {
        item.choices.forEach(choice => {
          all_nodes[+choice.id] = choice
          max_id = Math.max(max_id, +choice.id)
        })
      }
    })
  return { all_nodes, workflow, max_id }
}

export default function ChatThread() {
  const { data: session } = Auth.useSession()
  const router = useRouter()
  const [state, setState] = React.useState({} as Record<number, string>)
  const [message, setMessage] = React.useState('')
  const { data: messages, mutate } = useAPIQuery(GPTAssistantThread2MessagesList, { thread_id: router.query.thread_id as string })
  const { trigger, isMutating } = useAPIMutation(GPTAssistantThread2Message, { thread_id: router.query.thread_id as string })
  const { trigger: triggerDelete } = useAPIMutation(GPTAssistantThread2Delete, { thread_id: router.query.thread_id as string })
  const submit = React.useCallback(async (body: { message: string } | { step: { id: number, value?: string } }) => {
    const newMessages = await trigger({ body })
    await mutate(messages => [...messages ?? [], ...newMessages?.data ?? []])
    if ('step' in body && body.step.value) {
      setState({ [body.step.id]: body.step.value })
    }
    if (newMessages) {
      const lastMessage = newMessages.data[newMessages.data.length-1]
      if (lastMessage.role === 'assistant') {
        const lastMessageContent = lastMessage.content[lastMessage.content.length-1]
        if (lastMessageContent.type === 'text') {
          const value = TryJSONParse(lastMessageContent.text.value)
          if (value.suggestions.length === 1) {
            await submit({ step: value.suggestions[0] })
          }
        }
      }
    }
  }, [trigger])
  const playbookState = React.useMemo(() => messages ? AssembleState(
    messages
      .filter(msg => msg.role === 'user')
      .flatMap(msg => msg.content)
      .map(msgContent => msgContent.type === 'text' ? msgContent.text.value : undefined)
  ) : undefined, [messages])
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
                onClick={evt => {
                  submit({ message: suggestion })
                }}
              >{suggestion}</button>
            )
          })}
        </div>
        {messages?.map((message, i) =>
          <React.Fragment key={i}>
            {message.content.map((content, j) => {
              if (content.type === 'image_file') return null
              const value = TryJSONParse(content.text.value)
              const component = (typeof value === 'object' && 'step' in value) ? playbookState?.all_nodes[value.step.id] : undefined
              return (
                <React.Fragment key={j}>
                  {/*content.text.value*/}
                  {typeof value === 'string' ? <Message role={message.role} session={session}>{value}</Message>
                    : typeof value === 'object' && 'message' in value ? <Message role={message.role} session={session}>{value.message}</Message>
                    : null}
                  {component ? <Component component={component} state={state} setState={setState} /> : null}
                  {typeof value === 'object' && 'suggestions' in value && value.suggestions.length > 1 ?
                    <div className="flex flex-row flex-wrap justify-center gap-2 place-self-center">
                      {value.suggestions.map((suggestion: any) => {
                        const suggestionNode = playbookState?.all_nodes[suggestion.id] ? krg.getProcessNode(playbookState.all_nodes[suggestion.id].name) : undefined
                        return (
                          <div key={suggestion.id} className="tooltip" data-tip={suggestionNode?.meta?.description}>
                            <button
                              className="btn btn-ghost border border-primary btn-rounded rounded-lg btn-sm"
                              onClick={evt => {submit({ step: suggestion })}}
                            >
                              {suggestionNode ? <>
                              {suggestionNode.meta.label} {suggestion.likelihood !== undefined ? `(${suggestion.likelihood.toPrecision(2)})` : null}
                              </> : JSON.stringify(suggestion)}
                            </button>
                          </div>
                        )
                      })}
                    </div>
                    : null}
                </React.Fragment>
              )
            })}
          </React.Fragment>
        )}
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
