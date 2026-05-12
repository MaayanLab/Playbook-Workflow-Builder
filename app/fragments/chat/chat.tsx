import React from 'react'
import dynamic from 'next/dynamic'
import { GPTAssistantCreate, GPTAssistantMessage, GPTAssistantMessagesList } from "@/app/api/client"
import classNames from 'classnames'
import usePublicUrl from '@/utils/next-public-url'
import { z } from 'zod'

import { useAPIMutation, useAPIQuery } from "@/core/api/client"
import * as Auth from 'next-auth/react'

import krg from '@/app/krg'
import { useExRouter } from '@/app/fragments/ex-router'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { useFPL } from '../metapath'
import { StoryProvider } from '../story'
import { Waypoint, useWaypoints } from '@/app/components/waypoint'
import { Breadcrumbs } from '../breadcrumbs'
import { DataBreadcrumb, ProcessBreadcrumb } from '@/app/fragments/graph/breadcrumb'
import { close_icon, extend_icon, fullscreen_icon, func_icon, start_icon, variable_icon } from '@/icons'
import ReportButton from '../graph/report-button'
import { UnauthorizedError } from '@/spec/error'
import Link from 'next/link'

const Icon = dynamic(() => import('@/app/components/icon'))
const Cell = dynamic(() => import('@/app/fragments/report/cell'))
const Message = dynamic(() => import('@/app/fragments/chat/message'))
const SessionStatus = dynamic(() => import('@/app/fragments/session-status'))

export default function Page({ mode, session_id, graph_id, node_id, embedded = false }: { mode: string, session_id?: string, graph_id?: string, node_id?: string, embedded?: boolean }) {
  const router = useExRouter()
  const publicUrl = usePublicUrl()
  const [message, setMessage] = React.useState('')
  const { data: session } = Auth.useSession()
  const [thread_id, setThread_id] = React.useState<string | undefined>()
  const { data: { messages, fpl } = { messages: undefined, fpl: null }, mutate } = useAPIQuery(GPTAssistantMessagesList, () => thread_id ? { thread_id } : null, {
    refreshInterval: 1000,
    onSuccess(data) {
      if (data?.fpl && fpl !== data.fpl) {
        router.push(`${session_id ? `/session/${session_id}` : ''}/${mode}${mode !== 'chat' ? `/${data.fpl}` : ''}?thread_id=${thread_id}&message=`, undefined, { shallow: true, scroll: false })
      }
    }
  })
  const createMessage = useAPIMutation(GPTAssistantMessage, { thread_id })
  const createChat = useAPIMutation(GPTAssistantCreate)
  const submit = React.useCallback(async ({ thread_id, ...body }: { message: string, graph_id?: string, node_id?: string, thread_id?: string }) => {
    try {
      let thread_id_: string
      if (!thread_id) {
        thread_id_ = z.string().parse(await createChat.trigger())
        setThread_id(thread_id_)
      } else {
        thread_id_ = thread_id
      }
      const query = { thread_id: thread_id_ }
      await mutate((current) => ({ messages: [
        ...current?.messages ?? [],
        { id: '', role: 'user', content: body.message, fpl: graph_id ?? null, created: new Date(), thread: thread_id_, feedback: null },
      ], fpl: graph_id ?? null }), { revalidate: false })
      const res = await createMessage.trigger({ query, body })
      await mutate((current) => ({
        messages: array.unique([
          ...(current?.messages ?? []).slice(0, -1),
          ...res?.messages ?? [],
        ], x => x.id),
        fpl: res?.fpl ?? null,
      }), { revalidate: false })
    } catch (e: any) {
      if (UnauthorizedError.isinstance(e)) {
        const qs = new URLSearchParams(window.location.search)
        qs.set('message', body.message)
        Auth.signIn(undefined, { callbackUrl: `${window.location.pathname}?${qs.toString()}` })
      }
    }
  }, [])
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
  React.useEffect(() => {
    setThread_id(router.query.thread_id as string | undefined)
  }, [router.query.thread_id])
  React.useEffect(() => {
    if (typeof router.query.message === 'string') {
      submit({ message: router.query.message, graph_id, node_id })
    }
  }, [router.query.message])
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
            <ReportButton session_id={session_id} graph_id={fpl ?? 'start'} />
            </Waypoint>
            </>
          : null}
          <div className={classNames("m-2 grow flex flex-col justify-between overflow-hidden")}>
            {/* reversing elemens here is intentional since it keeps the scroll bar at the bottom */}
            <div className={classNames('flex flex-col-reverse px-2 bg-white dark:bg-current overflow-auto')}>
              {createMessage.error ?
                <Message role="assistant" session={session}>
                  I experienced an error: {createMessage.error instanceof Error ? createMessage.error.message : 'unknown error'}
                </Message>
                : null}
              {createMessage.isMutating ?
                <Message role="assistant" session={session}>
                  <span className="loading loading-dots loading-lg mt-2"></span>
                </Message>
                : null}
              {messages?.toReversed().map((message, i) => {
                const head = !embedded && 'fpl' in message && message.fpl && fpl_to_metapath[message.fpl]
                return (
                  <React.Fragment key={i}>
                    {head ?
                      <Cell
                        session_id={session_id}
                        krg={krg}
                        id={fpl ?? ''}
                        head={head}
                        cellMetadata={{ [head.id]: head.cell_metadata ?? { id: head.id, label: '', description: '', data_visible: true, process_visible: true } }}
                        setCellMetadata={() => {}}
                      />
                    : null}
                    <Message
                      thread_id={thread_id}
                      message_id={message.id}
                      role={message.role}
                      session={session}
                      embedded={embedded}
                    >
                      {message.content}
                    </Message>
                  </React.Fragment>
                )
              })}
              {(!graph_id || graph_id === 'start') && <div className="flex flex-row flex-wrap justify-center gap-2 place-self-center">
                {[
                  'Show me the expression of ACE2 in healthy human tissues from GTEx',
                  'Find drugs from the LINCS L1000 Chemical Perturbations that up regulate STAT3',
                  `I'd like to analyze the GEO study GSE301503`,
                ].map((suggestion, i) => {
                  return (
                    <button
                      key={i}
                      className="btn btn-primary btn-sm border border-primary btn-rounded rounded-lg"
                      onClick={evt => {submit({ message: suggestion, graph_id, node_id, thread_id })}}
                    >Example {i+1}</button>
                  )
                })}
              </div>}
              <Message role="welcome" session={session}>
                Hi! What would you like to build? Describe your idea and I'll generate the workflow for you.
              </Message>
              <div className={classNames("flex-grow max-w-none flex flex-col justify-center items-center")}>
                <img
                  className="w-32"
                  src={`${publicUrl}/PWB-logo.svg`}
                />
                <div className="prose text-center">
                  <h4>Text to Workflow</h4>
                </div>
              </div>
            </div>
            <Message role="user" session={session}>
              <form
                className="flex flex-col"
                onSubmit={async (evt) => {
                  evt.preventDefault()
                  const currentMessage = message
                  setMessage(() => '')
                  await submit({ message: currentMessage, graph_id, node_id, thread_id })
                }}
              >
                <textarea
                  className="textarea w-full bg-white"
                  rows={3}
                  placeholder="Type your questions here"
                  value={message}
                  onChange={evt => setMessage(() => evt.target.value)}
                  onKeyDown={async (evt) => {
                    if (evt.shiftKey && evt.key === 'Enter') {
                      evt.preventDefault()
                      evt.stopPropagation()
                      const currentMessage = message
                      setMessage(() => '')
                      await submit({ message: currentMessage, graph_id, node_id, thread_id })
                    }
                  }}
                />
                <button type="submit" className="btn btn-sm place-self-end" disabled={!message}>Send</button>
              </form>
            </Message>
          </div>
        </StoryProvider>
      </SessionStatus>
    </>
  )
}
