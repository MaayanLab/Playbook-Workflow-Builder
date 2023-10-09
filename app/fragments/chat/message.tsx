import React from 'react'
import dynamic from 'next/dynamic'

import classNames from 'classnames'
import { ReactMarkdown } from 'react-markdown/lib/react-markdown'
import krg from '@/app/krg'
import * as dict from '@/utils/dict'
import type { Session } from 'next-auth'

const UserAvatar = dynamic(() => import('@/app/fragments/playbook/avatar'))

type Component = {
  id: number,
  inputs: Record<string, number>,
  type: string,
  description: string,
  data?: any,
  likelihood?: number,
}

function Component({ state, setState, component }: {
  state: Record<number, any>, setState: React.Dispatch<React.SetStateAction<Record<number, any>>>,
  component: Component,
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
      formData.append(arg, JSON.stringify(state[component.inputs[arg]]))
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

export default function Message({
  session, role, children,
  submit, state, setState,
}: React.PropsWithChildren<{
  session: Session | null, role: string,
  submit: (message: { role: string, content: string }) => void,
  state: Record<number, string>, setState: React.Dispatch<React.SetStateAction<Record<number, string>>>,
}>) {
  const { component, suggestions }: { component?: Component, suggestions?: Component[] } = React.useMemo(() => {
    if (typeof children === 'string') {
      const m = /```(component|suggestions)\n(.+)\n```/g.exec(children)
      if (m) {
        return { [m[1]]: JSON.parse(m[2]) }
      }
    }
    return {}
  }, [children])
  if (suggestions) {
    return (
      <div className="flex flex-row gap-2 place-self-center">
      {suggestions.map(suggestion => {
        const suggestionNode = krg.getProcessNode(suggestion.type)
        return (
          <div key={suggestion.id} className="tooltip" data-tip={suggestionNode.meta.description}>
            <button
              className="btn btn-ghost border border-primary btn-rounded rounded-lg btn-sm"
              onClick={evt => {
                submit({ role: 'assistant', content: '```component\n' + JSON.stringify(suggestion) + '\n```' })
              }}
            >
              {suggestionNode.meta.label} {suggestion.likelihood !== undefined ? `(${suggestion.likelihood.toPrecision(2)})` : null}
            </button>
          </div>
        )
      })}
    </div>
    )
  }
  return (
    <>
      <div className={classNames('chat', { 'chat-end': role === 'user', 'chat-start': role !== 'user', 'hidden': role === 'system' })}>
        <div className="chat-image btn btn-ghost btn-circle avatar placeholder">
          <div className="bg-neutral-focus text-neutral-content rounded-full w-16">
            {role !== 'user' ? "PWB" : <UserAvatar session={session} />}
          </div>
        </div>
        <div className={classNames('chat-bubble rounded-xl w-full prose', { 'chat-bubble-primary': role === 'user', 'chat-bubble-secondary': role === 'assistant' || role === 'welcome', 'chat-bubble-error': role === 'error' })}>
          {component ? <Component component={component} state={state} setState={setState} />
            : typeof children === 'string' ? <ReactMarkdown>{children}</ReactMarkdown>
            : children}
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
