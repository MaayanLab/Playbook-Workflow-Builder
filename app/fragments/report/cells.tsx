import React from 'react'
import dynamic from 'next/dynamic'
import type KRG from '@/core/KRG'
import { type Metapath, useFPL } from '@/app/fragments/metapath'
import { StoryProvider } from '@/app/fragments/story'
import { useAPIMutation, useAPIQuery } from '@/core/api/client'
import { UpdateUserPlaybook, DeleteUserPlaybook, PublishUserPlaybook, UserPlaybook } from '@/app/api/client'
import * as dict from '@/utils/dict'
import { useRouter } from 'next/router'
import { Breadcrumbs } from '../breadcrumbs'
import { DataBreadcrumb, ProcessBreadcrumb } from '@/app/fragments/graph/breadcrumb'
import { extend_icon, func_icon, start_icon, variable_icon } from '@/icons'
import { Waypoint, useWaypoints } from '@/app/components/waypoint'

const Introduction = dynamic(() => import('@/app/fragments/report/introduction'))
const Cell = dynamic(() => import('@/app/fragments/report/cell'))
const SessionStatus = dynamic(() => import('@/app/fragments/session-status'))
const ImportButton = dynamic(() => import('@/app/fragments/graph/import-button'))
const CAVATICAButton = dynamic(() => import('@/app/fragments/graph/cavatica-button'))
const RestartButton = dynamic(() => import('@/app/fragments/graph/restart-button'))
const GraphButton = dynamic(() => import('@/app/fragments/report/graph-button'))
const Chat = dynamic(() => import('@/app/fragments/chat/chat'))

export default function Cells({ session_id, thread, krg, id }: { session_id?: string, thread?: string, krg: KRG, id: string }) {
  const router = useRouter()
  const { data: metapath } = useFPL(id)
  const { data: userPlaybook, mutate: mutateUserPlaybook } = useAPIQuery(UserPlaybook, { id })
  const data = React.useMemo(() => metapath ? ({ metapath, userPlaybook }) : undefined, [metapath, userPlaybook])
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
  }, [data])
  const process_to_step = React.useMemo(() => metapath ? dict.init(metapath.map(h => ({ key: h.process.id, value: `${h.id}:${h.process.id}` }))) : {}, [metapath])
  const head = React.useMemo(() => metapath ? metapath[metapath.length - 1] : undefined, [metapath])
  const { waypoints, scrollTo } = useWaypoints()
  if (!data || !playbookMetadata || !metapath) return null
  return (
    <div className="flex flex-col py-4 gap-2">
      {thread ? <Chat key={id} thread_id={thread} embedded /> : null}
      <SessionStatus session_id={session_id}>
        <StoryProvider krg={krg} metapath={data.metapath}>
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
                index={metapath.length * 2 + 1}
                id="extend"
                label="Extend"
                active={false}
                icon={extend_icon}
                parents={[head ? `${head.id}:${head.process.id}` : `start`]}
                onClick={() => {
                  router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${id}/extend`)
                }}
              />
            </Breadcrumbs>
            <ImportButton session_id={session_id} />
            <CAVATICAButton session_id={session_id} />
            <RestartButton session_id={session_id} />
            <GraphButton session_id={session_id} graph_id={id} />
          </Waypoint>
          <Waypoint id="start">
            <Introduction
              session_id={session_id}
              id={id}
              error={null}
              krg={krg}
              metapath={metapath}
              playbookMetadata={playbookMetadata}
              userPlaybook={userPlaybook}
              setPlaybookMetadata={setPlaybookMetadata}
              updateRequired={updateRequired}
              toggleSave={() => {
                if (updateRequired || !userPlaybook) {
                  const { id: playbookMetadataId, ...playbook_metadata } = playbookMetadata
                  const cell_metadata = dict.init(dict.items(cellMetadata).map(({ key, value }) => {
                    const { id, ...meta } = value
                    return { key, value: id ? { id } : meta }
                  }))
                  updateUserPlaybook({
                    query: { id },
                    body: {
                      user_playbook: { public: userPlaybook?.public || false },
                      playbook_metadata: playbookMetadataId ? { id: playbookMetadataId } : playbook_metadata,
                      cell_metadata,
                    },
                  }).then(id => {
                    mutateUserPlaybook({ public: userPlaybook?.public || false })
                    setPlaybookMetadata(metadata => ({ ...metadata, id: data.metapath[data.metapath.length-1].playbook_metadata?.id || '' }))
                    setUpdateRequired(false)
                    router.push(`${session_id ? `/session/${session_id}` : ''}/report/${id}`, undefined, { shallow: true, scroll: false })
                  })
                } else {
                  deleteUserPlaybook({
                    query: { id },
                    body: {},
                  }).then(id => {
                    mutateUserPlaybook(undefined)
                    setUpdateRequired(false)
                  })
                }
              }}
              togglePublic={() => {
                if (!updateRequired && userPlaybook) {
                  const publicPlaybook = !userPlaybook.public
                  publishUserPlaybook({
                    query: { id },
                    body: { public: publicPlaybook },
                  }).then(id => {
                    mutateUserPlaybook({ public: publicPlaybook })
                    setUpdateRequired(false)
                  })
                }
              }}
            />
          </Waypoint>
          {(data.metapath||[]).filter(head => head.id in cellMetadata).map(head => (
            <Cell
              key={`${head.id}-${cellMetadata[head.id]?.process_visible}-${cellMetadata[head.id]?.data_visible}`}
              session_id={session_id}
              krg={krg}
              id={id}
              head={head}
              cellMetadata={cellMetadata}
              setCellMetadata={setCellMetadata}
            />
          ))}
        </StoryProvider>
      </SessionStatus>
    </div>
  )
}