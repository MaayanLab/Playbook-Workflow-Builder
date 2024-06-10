import React from 'react'
import dynamic from 'next/dynamic'
import { GPTAssistantMessage, GPTAssistantMessagesList } from "@/app/api/client"
import classNames from 'classnames'
import usePublicUrl from '@/utils/next-public-url'

import { useAPIMutation, useAPIQuery } from "@/core/api/client"
import * as Auth from 'next-auth/react'

import krg from '@/app/krg'
import { AssembleState } from '@/app/api/v1/chat/utils'
import { useRouter } from 'next/router'
import * as dict from '@/utils/dict'
import { useFPL } from '../metapath'
import { StoryProvider } from '../story'
import { Waypoint, useWaypoints } from '@/app/components/waypoint'
import { Breadcrumbs } from '../breadcrumbs'
import { DataBreadcrumb, ProcessBreadcrumb } from '@/app/fragments/graph/breadcrumb'
import { close_icon, extend_icon, fullscreen_icon, func_icon, maximize_icon, minimize_icon, start_icon, variable_icon } from '@/icons'
import ReportButton from '../graph/report-button'
import Link from 'next/link'

const Cell = dynamic(() => import('@/app/fragments/report/cell'))
const Message = dynamic(() => import('@/app/fragments/chat/message'))
const Icon = dynamic(() => import('@/app/components/icon'))
const SessionStatus = dynamic(() => import('@/app/fragments/session-status'))

export default function Page({ thread_id, session_id, embedded = false }: { thread_id: string, session_id?: string, embedded?: boolean }) {
  const router = useRouter()
  const publicUrl = usePublicUrl()
  const [message, setMessage] = React.useState('')
  const [collapse, setCollapse] = React.useState(false)
  const { data: session } = Auth.useSession({ required: true })
  const { data: { messages, fpl } = { messages: undefined, fpl: null }, mutate } = useAPIQuery(GPTAssistantMessagesList, { thread_id })
  const { trigger, isMutating } = useAPIMutation(GPTAssistantMessage, { thread_id })
  // const { trigger: triggerDelete } = useAPIMutation(GPTAssistantDelete, { thread_id })
  const playbookState = React.useMemo(() => messages ? AssembleState(messages, { with_value: true }) : undefined, [messages])
  const submit = React.useCallback(async (body: { message: string } | { step: { id: number, value?: string } }) => {
    const res = await trigger({ body })
    if (res?.fpl) {
      router.push(`${session_id ? `/session/${session_id}` : ''}/report/${res.fpl}?thread=${thread_id}`)
    } else {
      await mutate((current) => ({ messages: [...current?.messages ?? [], ...res?.messages ?? []], fpl: null }))
    }
  }, [trigger, session_id, thread_id])
  const { data: metapath } = useFPL(fpl ? fpl : undefined)
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
      <style jsx>{`
        /* From the newer version of DaisyUI */
        .loading {
          pointer-events: none;
          display: inline-block;
          aspect-ratio: 1/1;
          width: 1.5rem;
          background-color: currentColor;
          -webkit-mask-size: 100%;
          mask-size: 100%;
          -webkit-mask-repeat: no-repeat;
          mask-repeat: no-repeat;
          -webkit-mask-position: center;
          mask-position: center;
          -webkit-mask-image: url(data:image/svg+xml;base64,PHN2ZyB3aWR0aD0nMjQnIGhlaWdodD0nMjQnIHN0cm9rZT0nIzAwMCcgdmlld0JveD0nMCAwIDI0IDI0JyB4bWxucz0naHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmcnPjxzdHlsZT4uc3Bpbm5lcl9WOG0xe3RyYW5zZm9ybS1vcmlnaW46Y2VudGVyO2FuaW1hdGlvbjpzcGlubmVyX3pLb2EgMnMgbGluZWFyIGluZmluaXRlfS5zcGlubmVyX1Y4bTEgY2lyY2xle3N0cm9rZS1saW5lY2FwOnJvdW5kO2FuaW1hdGlvbjpzcGlubmVyX1lwWlMgMS41cyBlYXNlLW91dCBpbmZpbml0ZX1Aa2V5ZnJhbWVzIHNwaW5uZXJfektvYXsxMDAle3RyYW5zZm9ybTpyb3RhdGUoMzYwZGVnKX19QGtleWZyYW1lcyBzcGlubmVyX1lwWlN7MCV7c3Ryb2tlLWRhc2hhcnJheTowIDE1MDtzdHJva2UtZGFzaG9mZnNldDowfTQ3LjUle3N0cm9rZS1kYXNoYXJyYXk6NDIgMTUwO3N0cm9rZS1kYXNob2Zmc2V0Oi0xNn05NSUsMTAwJXtzdHJva2UtZGFzaGFycmF5OjQyIDE1MDtzdHJva2UtZGFzaG9mZnNldDotNTl9fTwvc3R5bGU+PGcgY2xhc3M9J3NwaW5uZXJfVjhtMSc+PGNpcmNsZSBjeD0nMTInIGN5PScxMicgcj0nOS41JyBmaWxsPSdub25lJyBzdHJva2Utd2lkdGg9JzMnPjwvY2lyY2xlPjwvZz48L3N2Zz4=);
          mask-image: url(data:image/svg+xml;base64,PHN2ZyB3aWR0aD0nMjQnIGhlaWdodD0nMjQnIHN0cm9rZT0nIzAwMCcgdmlld0JveD0nMCAwIDI0IDI0JyB4bWxucz0naHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmcnPjxzdHlsZT4uc3Bpbm5lcl9WOG0xe3RyYW5zZm9ybS1vcmlnaW46Y2VudGVyO2FuaW1hdGlvbjpzcGlubmVyX3pLb2EgMnMgbGluZWFyIGluZmluaXRlfS5zcGlubmVyX1Y4bTEgY2lyY2xle3N0cm9rZS1saW5lY2FwOnJvdW5kO2FuaW1hdGlvbjpzcGlubmVyX1lwWlMgMS41cyBlYXNlLW91dCBpbmZpbml0ZX1Aa2V5ZnJhbWVzIHNwaW5uZXJfektvYXsxMDAle3RyYW5zZm9ybTpyb3RhdGUoMzYwZGVnKX19QGtleWZyYW1lcyBzcGlubmVyX1lwWlN7MCV7c3Ryb2tlLWRhc2hhcnJheTowIDE1MDtzdHJva2UtZGFzaG9mZnNldDowfTQ3LjUle3N0cm9rZS1kYXNoYXJyYXk6NDIgMTUwO3N0cm9rZS1kYXNob2Zmc2V0Oi0xNn05NSUsMTAwJXtzdHJva2UtZGFzaGFycmF5OjQyIDE1MDtzdHJva2UtZGFzaG9mZnNldDotNTl9fTwvc3R5bGU+PGcgY2xhc3M9J3NwaW5uZXJfVjhtMSc+PGNpcmNsZSBjeD0nMTInIGN5PScxMicgcj0nOS41JyBmaWxsPSdub25lJyBzdHJva2Utd2lkdGg9JzMnPjwvY2lyY2xlPjwvZz48L3N2Zz4=)
        }
        .loading-dots {
          -webkit-mask-image: url(data:image/svg+xml;base64,PHN2ZyB3aWR0aD0nMjQnIGhlaWdodD0nMjQnIHZpZXdCb3g9JzAgMCAyNCAyNCcgeG1sbnM9J2h0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnJz48c3R5bGU+LnNwaW5uZXJfcU04M3thbmltYXRpb246c3Bpbm5lcl84SFFHIDEuMDVzIGluZmluaXRlfS5zcGlubmVyX29YUHJ7YW5pbWF0aW9uLWRlbGF5Oi4xc30uc3Bpbm5lcl9aVExme2FuaW1hdGlvbi1kZWxheTouMnN9QGtleWZyYW1lcyBzcGlubmVyXzhIUUd7MCUsNTcuMTQle2FuaW1hdGlvbi10aW1pbmctZnVuY3Rpb246Y3ViaWMtYmV6aWVyKDAuMzMsLjY2LC42NiwxKTt0cmFuc2Zvcm06dHJhbnNsYXRlKDApfTI4LjU3JXthbmltYXRpb24tdGltaW5nLWZ1bmN0aW9uOmN1YmljLWJlemllcigwLjMzLDAsLjY2LC4zMyk7dHJhbnNmb3JtOnRyYW5zbGF0ZVkoLTZweCl9MTAwJXt0cmFuc2Zvcm06dHJhbnNsYXRlKDApfX08L3N0eWxlPjxjaXJjbGUgY2xhc3M9J3NwaW5uZXJfcU04MycgY3g9JzQnIGN5PScxMicgcj0nMycvPjxjaXJjbGUgY2xhc3M9J3NwaW5uZXJfcU04MyBzcGlubmVyX29YUHInIGN4PScxMicgY3k9JzEyJyByPSczJy8+PGNpcmNsZSBjbGFzcz0nc3Bpbm5lcl9xTTgzIHNwaW5uZXJfWlRMZicgY3g9JzIwJyBjeT0nMTInIHI9JzMnLz48L3N2Zz4=);
          mask-image: url(data:image/svg+xml;base64,PHN2ZyB3aWR0aD0nMjQnIGhlaWdodD0nMjQnIHZpZXdCb3g9JzAgMCAyNCAyNCcgeG1sbnM9J2h0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnJz48c3R5bGU+LnNwaW5uZXJfcU04M3thbmltYXRpb246c3Bpbm5lcl84SFFHIDEuMDVzIGluZmluaXRlfS5zcGlubmVyX29YUHJ7YW5pbWF0aW9uLWRlbGF5Oi4xc30uc3Bpbm5lcl9aVExme2FuaW1hdGlvbi1kZWxheTouMnN9QGtleWZyYW1lcyBzcGlubmVyXzhIUUd7MCUsNTcuMTQle2FuaW1hdGlvbi10aW1pbmctZnVuY3Rpb246Y3ViaWMtYmV6aWVyKDAuMzMsLjY2LC42NiwxKTt0cmFuc2Zvcm06dHJhbnNsYXRlKDApfTI4LjU3JXthbmltYXRpb24tdGltaW5nLWZ1bmN0aW9uOmN1YmljLWJlemllcigwLjMzLDAsLjY2LC4zMyk7dHJhbnNmb3JtOnRyYW5zbGF0ZVkoLTZweCl9MTAwJXt0cmFuc2Zvcm06dHJhbnNsYXRlKDApfX08L3N0eWxlPjxjaXJjbGUgY2xhc3M9J3NwaW5uZXJfcU04MycgY3g9JzQnIGN5PScxMicgcj0nMycvPjxjaXJjbGUgY2xhc3M9J3NwaW5uZXJfcU04MyBzcGlubmVyX29YUHInIGN4PScxMicgY3k9JzEyJyByPSczJy8+PGNpcmNsZSBjbGFzcz0nc3Bpbm5lcl9xTTgzIHNwaW5uZXJfWlRMZicgY3g9JzIwJyBjeT0nMTInIHI9JzMnLz48L3N2Zz4=)
        }
        .loading-lg {
          width: 2.5rem;
        }
      `}</style>
      <SessionStatus session_id={session_id}>
        <StoryProvider krg={krg} metapath={metapath ?? []}>
          {!embedded && metapath ?
          <>
            <Waypoint id="head" className="sticky top-0 left-0 z-20 bg-white dark:bg-current w-full flex flex-row place-items-center">
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
              {metapath.flatMap((step, i) => {
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
                      // setCellMetadata((cellMetadata) => ({ ...cellMetadata, [step.id]: { ...cellMetadata[step.id], process_visible: true, id: '' } }))
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
                      // setCellMetadata((cellMetadata) => ({ ...cellMetadata, [step.id]: { ...cellMetadata[step.id], data_visible: true, id: '' } }))
                      scrollTo(`${step.id}:data`)
                    }}
                  />,
                ]
              })}
              <ProcessBreadcrumb
                key="extend"
                index={metapath.length * 2 + 1}
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
            <ReportButton session_id={session_id} graph_id={fpl ?? 'start'} thread_id={thread_id} />
            </Waypoint>
            </>
            : null}
          <div className={classNames('flex flex-col', {"absolute top-0 left-1/2 w-1/2 z-30 h-screen max-h-screen mb-5 mr-5 pr-10 bg-transparent justify-end overflow-hidden pointer-events-none": embedded})}>
            <div className={classNames('flex flex-col px-2 bg-white dark:bg-current', { 'border rounded-xl border-black mt-48 pointer-events-auto overflow-hidden': embedded})}>
              <div className={classNames('flex-grow flex flex-row my-1 gap-2 bg-white dark:bg-current p-2', {'hidden': !embedded})}>
                <div className="prose"><h3>Text to Workflow</h3></div>
                <div className="flex-grow">&nbsp;</div>
                <button onClick={() => {setCollapse(c => !c)}}>
                  <Icon icon={collapse ? maximize_icon : minimize_icon} className="fill-black dark:fill-white" />
                </button>
                {!collapse ? <>
                  <Link href={`/chat/${thread_id}`}>
                    <button>
                      <Icon icon={fullscreen_icon} className="fill-black dark:fill-white" />
                    </button>
                  </Link>
                  <Link href={`/report/${fpl}`} shallow>
                    <button className="bg-red-500">
                      <Icon icon={close_icon} className="fill-black dark:fill-white" />
                    </button>
                  </Link>
                </> : null}
              </div>
              <div className={classNames('flex-grow flex flex-col overflow-hidden', {'bg-yellow-50': embedded, 'hidden': collapse })}>
                <div className='flex-grow flex flex-col-reverse overflow-y-auto overflow-x-hidden'>
                  <div>
                    <div className={classNames("flex-grow max-w-none flex flex-col justify-center items-center")}>
                      <img
                        className="w-32"
                        src={`${publicUrl}/PWB-logo.svg`}
                      />
                      <div className="prose"><h4>Playbook Workflow Builder Text to Workflow</h4></div>
                    </div>
                    <Message role="welcome" session={session}>
                      How can I help you today?
                    </Message>
                    <div className="flex flex-row flex-wrap justify-center gap-2 place-self-center">
                      {[
                        'Show me the expression of ACE2 in healthy human tissues from GTEx',
                        'Find drugs from the LINCS L1000 Chemical Perturbations that up regulate STAT3',
                      ].map((suggestion, i) => {
                        return (
                          <button
                            key={i}
                            className="btn btn-ghost border border-primary btn-rounded rounded-lg btn-sm bg-white"
                            onClick={evt => {submit({ message: suggestion })}}
                          >{suggestion}</button>
                        )
                      })}
                    </div>
                    {messages?.map((message, i) => {
                      const head = !embedded && 'fpl' in message && message.fpl && fpl_to_metapath[message.fpl]
                      return (
                        <React.Fragment key={i}>
                          {head ?
                            <Cell
                              key={message.fpl}
                              session_id={session_id}
                              krg={krg}
                              id={fpl ?? ''}
                              head={head}
                              cellMetadata={{ [head.id]: head.cell_metadata ?? { id: head.id, label: '', description: '', data_visible: true, process_visible: true } }}
                              setCellMetadata={() => {}}
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
                                      className="btn btn-ghost border border-primary btn-rounded rounded-lg btn-sm bg-white"
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
                        <span className="loading loading-dots loading-lg mt-2"></span>
                      </Message>
                      : null}
                  </div>
                </div>
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
            </div>
          </div>
        </StoryProvider>
      </SessionStatus>
    </>
  )
}
