import React from 'react'
import dynamic from 'next/dynamic'
import { type Metapath, useFPL } from '@/app/fragments/metapath'
import { StoryProvider } from '@/app/fragments/story'
import { UpdateUserPlaybook, DeleteUserPlaybook, PublishUserPlaybook } from '@/app/api/client'
import { GPTAssistantMessage, GPTAssistantMessageFeedback, GPTAssistantMessagesList } from "@/app/api/client"
import * as dict from '@/utils/dict'
import { useRouter } from 'next/router'
import { Breadcrumbs } from '../breadcrumbs'
import { DataBreadcrumb, ProcessBreadcrumb } from '@/app/fragments/graph/breadcrumb'
import { extend_icon, func_icon, start_icon, variable_icon } from '@/icons'
import { Waypoint, useWaypoints } from '@/app/components/waypoint'
import classNames from 'classnames'
import usePublicUrl from '@/utils/next-public-url'

import { useAPIMutation, useAPIQuery } from "@/core/api/client"
import * as Auth from 'next-auth/react'

import krg from '@/app/krg'
import type { Session } from 'next-auth'
import { AssembleState } from '@/app/api/v1/chat/utils'
import type { ReactMarkdownOptions } from 'react-markdown/lib/react-markdown'

const Introduction = dynamic(() => import('@/app/fragments/report/introduction'))
const Cell = dynamic(() => import('@/app/fragments/report/cell'))
const SessionStatus = dynamic(() => import('@/app/fragments/session-status'))
const ImportButton = dynamic(() => import('@/app/fragments/graph/import-button'))
const CAVATICAButton = dynamic(() => import('@/app/fragments/graph/cavatica-button'))
const RestartButton = dynamic(() => import('@/app/fragments/graph/restart-button'))
const GraphButton = dynamic(() => import('@/app/fragments/report/graph-button'))
const UserAvatar = dynamic(() => import('@/app/fragments/playbook/avatar'))
const ReactMarkdown = dynamic(() => import('react-markdown/lib/react-markdown').then(({ ReactMarkdown }) => ReactMarkdown as ((props: ReactMarkdownOptions) => React.ReactNode)), { ssr: false })

function Message({ thread_id, message_id, session, role, children }: React.PropsWithChildren<{ thread_id?: string, message_id?: string, session: Session | null, role: string }>) {
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

export default function Page({ thread_id, session_id }: { thread_id: string, session_id?: string }) {
  const publicUrl = usePublicUrl()
  const [message, setMessage] = React.useState('')
  const { data: session } = Auth.useSession()
  const { data: { messages, fpl } = { messages: undefined, fpl: null }, mutate } = useAPIQuery(GPTAssistantMessagesList, { thread_id })
  const { trigger, isMutating } = useAPIMutation(GPTAssistantMessage, { thread_id })
  // const { trigger: triggerDelete } = useAPIMutation(GPTAssistantDelete, { thread_id })
  const playbookState = React.useMemo(() => messages ? AssembleState(messages, { with_value: true }) : undefined, [messages])
  const submit = React.useCallback(async (body: { message: string } | { step: { id: number, value?: string } }) => {
    const res = await trigger({ body })
    await mutate((current) => ({ messages: [...current?.messages ?? [], ...res?.messages ?? []], fpl: res?.fpl ?? null }))
  }, [trigger])
  const router = useRouter()
  const { data: metapath } = useFPL(fpl ? fpl : undefined)
  const data = React.useMemo(() => metapath?.length ? ({ metapath, userPlaybook: undefined }) : undefined, [metapath])
  const { trigger: updateUserPlaybook } = useAPIMutation(UpdateUserPlaybook, undefined, { base: session_id ? `/api/socket/${session_id}` : '', throwOnError: true })
  const { trigger: publishUserPlaybook } = useAPIMutation(PublishUserPlaybook, undefined, { base: session_id ? `/api/socket/${session_id}` : '', throwOnError: true })
  const { trigger: deleteUserPlaybook } = useAPIMutation(DeleteUserPlaybook, undefined, { base: session_id ? `/api/socket/${session_id}` : '', throwOnError: true })
  const [cellMetadata, setCellMetadata] = React.useState({} as Record<string, Exclude<Metapath['cell_metadata'], null>>)
  const [playbookMetadata, setPlaybookMetadata] = React.useState({
    id: '',
    title: '',
    description: '',
    summary: 'auto',
    gpt_summary: '',
  } as Exclude<Metapath['playbook_metadata'], null>)
  const [userPlaybook, setUserPlaybook] = React.useState(undefined as undefined | { public: boolean })
  const [updateRequired, setUpdateRequired] = React.useState(false)
  React.useEffect(() => {
    if (!data) return
    setUpdateRequired(
      playbookMetadata?.id !== data.metapath[data.metapath.length-1].playbook_metadata?.id
      || data.metapath.some(cell => cell.cell_metadata?.id !== cellMetadata[cell.id]?.id)
    )
  }, [data, cellMetadata, playbookMetadata, userPlaybook])
  React.useEffect(() => {
    if (!data) return
    const {
      id = '',
      title = '',
      description = '',
      summary = 'auto',
      gpt_summary = '',
    } = data.metapath[data.metapath.length-1].playbook_metadata ?? {}
    setPlaybookMetadata({ id, title, description, summary, gpt_summary })
    setCellMetadata(dict.init(data.metapath.map((element, index) => {
      const {
        id = '',
        label = '',
        description = '',
        process_visible = false,
        data_visible = index+1 !== data.metapath.length-1,
      } = element.cell_metadata ?? {}
      return { key: element.id, value: { id, label, description, process_visible, data_visible } }
    })))
    setUserPlaybook(data.userPlaybook)
  }, [data])
  const { fpl_to_metapath, process_to_step } = React.useMemo(() => metapath ? {
    fpl_to_metapath: dict.init(metapath.map(h => ({ key: h.id, value: h }))),
    process_to_step: dict.init(metapath.map(h => ({ key: h.process.id, value: `${h.id}:${h.process.id}` }))),
  } : {
    fpl_to_metapath: {},
    process_to_step: {},
  }, [metapath])
  const head = React.useMemo(() => metapath ? metapath[metapath.length - 1] : undefined, [metapath])
  const { waypoints, scrollTo } = useWaypoints()
  return (
    <>
      <SessionStatus session_id={session_id}>
        <StoryProvider krg={krg} metapath={metapath ?? []}>
        <Waypoint id="head" className="sticky top-0 left-0 z-50 bg-white dark:bg-current w-full flex flex-row place-items-center">
          <Breadcrumbs>
            <DataBreadcrumb
              key="start"
              index={0}
              id="start"
              label="Start"
              active={waypoints.get('start')?.active !== false}
              icon={[start_icon]}
              parents={[]}
              onClick={() => {
                scrollTo('top')
              }}
            />
            {metapath?.flatMap((step, i) => {
              const process = krg.getProcessNode(step.process.type)
              if (process === undefined) return []
              return [
                <ProcessBreadcrumb
                  key={step.id}
                  index={i * 2 + 1}
                  id={step.id}
                  label={process.meta.label}
                  head={step}
                  active={false}
                  icon={process.meta.icon || [func_icon]}
                  parents={dict.isEmpty(step.process.inputs) ? ['start'] : dict.values(step.process.inputs).map(({ id }) => process_to_step[id])}
                  onClick={() => {
                    setCellMetadata((cellMetadata) => ({ ...cellMetadata, [step.id]: { ...cellMetadata[step.id], process_visible: true, id: '' } }))
                    scrollTo(`${step.id}:process`)
                  }}
                />,
                <DataBreadcrumb
                  key={`${step.id}:${step.process.id}`}
                  index={i * 2 + 2}
                  id={`${step.id}:${step.process.id}`}
                  label={process.output.meta.label}
                  head={step}
                  active={!!waypoints.get(`${step.id}:data`)?.active}
                  icon={process.output.meta.icon || [variable_icon]}
                  parents={[step.id]}
                  onClick={() => {
                    setCellMetadata((cellMetadata) => ({ ...cellMetadata, [step.id]: { ...cellMetadata[step.id], data_visible: true, id: '' } }))
                    scrollTo(`${step.id}:data`)
                  }}
                />,
              ]
            })}
            <ProcessBreadcrumb
              key="extend"
              index={(metapath?.length ?? 0) * 2 + 1}
              id="extend"
              label="Extend"
              active={false}
              icon={extend_icon}
              parents={[head ? `${head.id}:${head.process.id}` : `start`]}
              onClick={() => {
                scrollTo(`bottom`)
              }}
            />
          </Breadcrumbs>
          <ImportButton session_id={session_id} />
          <CAVATICAButton session_id={session_id} />
          <RestartButton session_id={session_id} />
          {fpl ? <GraphButton session_id={session_id} graph_id={fpl} /> : null}
        </Waypoint>
        <div className="flex-grow flex flex-col gap-4">
          <div className={classNames("flex-grow max-w-none flex flex-col justify-center items-center", { 'hidden': isMutating || messages?.length })}>
            <img
              className="w-32"
              src={`${publicUrl}/PWB-logo.svg`}
            />
            <div className="prose"><h1>How can I help you today?</h1></div>
          </div>
          {data ?
            <Waypoint id="start">
              <Introduction
                session_id={session_id}
                id={fpl ?? 'start'}
                error={null}
                krg={krg}
                metapath={metapath ?? []}
                playbookMetadata={playbookMetadata}
                userPlaybook={userPlaybook}
                setPlaybookMetadata={setPlaybookMetadata}
                updateRequired={updateRequired}
                toggleSave={() => {
                  if (!fpl || !data) return
                  if (updateRequired || !userPlaybook) {
                    const { id: playbookMetadataId, ...playbook_metadata } = playbookMetadata
                    const cell_metadata = dict.init(dict.items(cellMetadata).map(({ key, value }) => {
                      const { id, ...meta } = value
                      return { key, value: id ? { id } : meta }
                    }))
                    updateUserPlaybook({
                      query: { id: fpl },
                      body: {
                        user_playbook: { public: userPlaybook?.public || false },
                        playbook_metadata: playbookMetadataId ? { id: playbookMetadataId } : playbook_metadata,
                        cell_metadata,
                      },
                    }).then(id => {
                      setUserPlaybook({ public: userPlaybook?.public || false })
                      setPlaybookMetadata(metadata => ({ ...metadata, id: data.metapath[data.metapath.length-1].playbook_metadata?.id || '' }))
                      setUpdateRequired(false)
                      router.push(`${session_id ? `/session/${session_id}` : ''}/report/${id}`, undefined, { shallow: true, scroll: false })
                    })
                  } else {
                    deleteUserPlaybook({
                      query: { id: fpl },
                      body: {},
                    }).then(id => {
                      setUserPlaybook(undefined)
                      setUpdateRequired(false)
                    })
                  }
                }}
                togglePublic={() => {
                  if (!fpl) return
                  if (!updateRequired && userPlaybook) {
                    const publicPlaybook = !userPlaybook.public
                    publishUserPlaybook({
                      query: { id: fpl },
                      body: { public: publicPlaybook },
                    }).then(id => {
                      setUserPlaybook({ public: publicPlaybook })
                      setUpdateRequired(false)
                    })
                  }
                }}
              />
            </Waypoint>
          : <>
            <Message role="welcome" session={session}>
              I'm an AI-powered chat assistant interface designed to help you access the functionality of the playbook workflow builder.
              Please start by asking your question of interest, and I'll try my best to help you answer it through the construction of a playbook workflow.
            </Message>
            <div className="flex flex-row flex-wrap justify-center gap-2 place-self-center prose">
              {[
                'Show me the expression of ACE2 in healthy human tissues from GTEx',
                'Find drugs from the LINCS L1000 Chemical Perturbations that up regulate STAT3',
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
          </>}
          {/* <div className="flex-grow prose max-w-none flex flex-row justify-between">
            <button
              type="button"
              className="btn btn-sm btn-error"
              onClick={() => {triggerDelete().finally(() => {router.push('/chat')})}}
            >Close</button>
          </div> */}
          {messages?.map((message, i) => {
            return (
              <React.Fragment key={i}>
                {'fpl' in message && message.fpl && fpl_to_metapath[message.fpl] ?
                  <Cell
                    key={`${message.fpl}-${cellMetadata[message.fpl]?.process_visible}-${cellMetadata[message.fpl]?.data_visible}`}
                    session_id={session_id}
                    krg={krg}
                    id={fpl ?? 'start'}
                    head={fpl_to_metapath[message.fpl]}
                    cellMetadata={cellMetadata}
                    setCellMetadata={setCellMetadata}
                  />
                  : null}
                {'message' in message ?
                  <Message
                    thread_id={thread_id}
                    message_id={message.id}
                    role={message.role}
                    session={session}
                  >{message.message}</Message>
                  : null}
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
        </div>
        </StoryProvider>
      </SessionStatus>
    </>
  )
}
