import React from 'react'
import dynamic from 'next/dynamic'
import classNames from 'classnames'
import usePublicUrl from '@/utils/next-public-url'
import { z } from 'zod'

import { useQuery, useMutation, useQueryClient, skipToken } from '@tanstack/react-query'
import { useTRPC, useTRPCClient, isTRPCClientError } from '@/app/trpc/client'
import * as Auth from '@/app/extensions/auth-client'

import { useFPL, useKRG } from '@/app/fragments/metapath'
import { useExRouter } from '@/app/fragments/ex-router'
import * as dict from '@/utils/dict'
import { StoryProvider } from '../story'
import { Waypoint, useWaypoints } from '@/app/components/waypoint'
import { Breadcrumbs } from '../breadcrumbs'
import { DataBreadcrumb, ProcessBreadcrumb } from '@/app/fragments/graph/breadcrumb'
import { extend_icon, func_icon, start_icon, variable_icon } from '@/icons'
import ReportButton from '../graph/report-button'

const Icon = dynamic(() => import('@/app/components/icon'))
const Cell = dynamic(() => import('@/app/fragments/report/cell'))
const Message = dynamic(() => import('@/app/fragments/chat/message'))
const SessionStatus = dynamic(() => import('@/app/fragments/session-status'))

export default function Page({ mode, session_id, graph_id, node_id, thread_id: initialThreadId, embedded = false }: { mode: string, session_id?: string, graph_id?: string, node_id?: string, thread_id?: string, embedded?: boolean }) {
  const krg = useKRG()
  const router = useExRouter()
  const publicUrl = usePublicUrl()
  const [message, setMessage] = React.useState('')
  const { data: session } = Auth.useSession()
  const [thread_id, setThread_id] = React.useState<string | undefined>(initialThreadId)
  const trpc = useTRPC()
  const trpcClient = useTRPCClient()
  const queryClient = useQueryClient()
  const { data, isFetched } = useQuery({
    ...trpc.chat.messages.queryOptions(thread_id ? { thread_id } : skipToken),
    refetchInterval: 1000,
  })
  const { messages, fpl } = data ?? { messages: undefined, fpl: null }
  // navigate to the workflow when the assistant moves the thread to a new one: the workflow the
  //  thread is on when it's first loaded is where it was, not a move (the person may be viewing
  //  any of the thread's workflows, or one in another mode, on purpose)
  const previousFpl = React.useRef<string | null | undefined>(undefined)
  // a different thread starts over
  React.useEffect(() => { previousFpl.current = undefined }, [thread_id])
  React.useEffect(() => {
    if (!isFetched) return
    if (previousFpl.current !== undefined && fpl && previousFpl.current !== fpl) {
      router.push(`${session_id ? `/session/${session_id}` : ''}/${mode}${mode !== 'chat' ? `/${fpl}` : ''}?thread_id=${thread_id}&message=`, undefined, { shallow: true, scroll: false })
    }
    previousFpl.current = fpl
  }, [isFetched, fpl])
  const createChat = useMutation(trpc.chat.create.mutationOptions())
  const [sendPending, setSendPending] = React.useState(false)
  const [sendError, setSendError] = React.useState<Error | null>(null)
  const submit = React.useCallback(async ({ thread_id, ...body }: { message: string, graph_id?: string, node_id?: string, thread_id?: string }) => {
    try {
      let thread_id_: string
      if (!thread_id) {
        thread_id_ = z.string().parse(await createChat.mutateAsync())
        setThread_id(thread_id_)
      } else {
        thread_id_ = thread_id
      }
      const queryKey = trpc.chat.messages.queryKey({ thread_id: thread_id_ })
      await queryClient.cancelQueries({ queryKey })
      queryClient.setQueryData(queryKey, (current) => ({ messages: [
        ...current?.messages ?? [],
        { id: '', role: 'user', content: body.message, fpl: graph_id ?? null, created: new Date().toISOString(), thread: thread_id_, feedback: null },
      ], fpl: graph_id ?? null }))
      setSendPending(true)
      setSendError(null)
      await new Promise<void>((resolve, reject) => {
        trpcClient.chat.send.subscribe({ thread_id: thread_id_, ...body }, {
          onData: (update) => {
            queryClient.setQueryData(queryKey, (current) => {
              const updatedIds = new Set(update.messages.map(message => message.id))
              return {
                messages: [
                  ...(current?.messages ?? []).filter(message => message.id !== '' && !updatedIds.has(message.id)),
                  ...update.messages,
                ],
                fpl: update.fpl ?? null,
              }
            })
          },
          onError: reject,
          onComplete: resolve,
        })
      })
    } catch (e: any) {
      if (isTRPCClientError(e, 'UNAUTHORIZED')) {
        const qs = new URLSearchParams(window.location.search)
        qs.set('message', body.message)
        Auth.signIn(undefined, { callbackUrl: `${window.location.pathname}?${qs.toString()}` })
      } else {
        setSendError(e instanceof Error ? e : new Error('unknown error'))
      }
    } finally {
      setSendPending(false)
    }
  }, [createChat, graph_id, queryClient, trpc, trpcClient])
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
    setThread_id((router.query.thread_id as string | undefined) ?? initialThreadId)
  }, [initialThreadId, router.query.thread_id])
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
              {sendError ?
                <Message role="assistant" session={session}>
                  I experienced an error: {sendError.message}
                </Message>
                : null}
              {sendPending ?
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
                      title={suggestion}
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
