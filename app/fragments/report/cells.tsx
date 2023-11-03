import React from 'react'
import dynamic from 'next/dynamic'
import type KRG from '@/core/KRG'
import { type Metapath, useFPL } from '@/app/fragments/metapath'
import { StoryProvider } from '@/app/fragments/story'
import { useAPIMutation } from '@/core/api/client'
import { UserPlaybook, UpdateUserPlaybook, DeleteUserPlaybook, PublishUserPlaybook } from '@/app/api/client'
import * as dict from '@/utils/dict'
import { useRouter } from 'next/router'
import { Breadcrumbs } from '../breadcrumbs'
import { DataBreadcrumb, ProcessBreadcrumb } from '../graph/breadcrumb'
import { extend_icon, func_icon, start_icon, variable_icon } from '@/icons'
import { Waypoint, useWaypoints } from '@/app/components/waypoint'

const Introduction = dynamic(() => import('@/app/fragments/report/introduction'))
const Cell = dynamic(() => import('@/app/fragments/report/cell'))
const SessionStatus = dynamic(() => import('@/app/fragments/session-status'))

export default function Cells({ session_id, krg, id }: { session_id?: string, krg: KRG, id: string }) {
  const router = useRouter()
  const { data: metapath } = useFPL(id)
  const data = React.useMemo(() => metapath ? ({ metapath, userPlaybook: undefined }) : undefined, [metapath])
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
  const process_to_step = React.useMemo(() => metapath ? dict.init(metapath.map(h => ({ key: h.process.id, value: `${h.id}:${h.process.id}` }))) : {}, [metapath])
  const head = React.useMemo(() => metapath ? metapath[metapath.length - 1] : undefined, [metapath])
  const waypoints = useWaypoints()
  if (!data || !playbookMetadata || !metapath) return null
  return (
    <div className="flex flex-col py-4 gap-2">
      <SessionStatus session_id={session_id}>
        <div className="sticky top-0 left-0 z-50 bg-white w-full">
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
                router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${id}/node/start`, undefined, { shallow: true })
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
                    router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${id}${id !== step.id ? `/node/${step.id}` : ''}`, undefined, { shallow: true })
                  }}
                />,
                <DataBreadcrumb
                  key={`${step.id}:${step.process.id}`}
                  index={i * 2 + 2}
                  id={`${step.id}:${step.process.id}`}
                  label={process.output.meta.label}
                  head={step}
                  active={!!waypoints.get(step.id)?.active}
                  icon={process.output.meta.icon || [variable_icon]}
                  parents={[step.id]}
                  onClick={() => {
                    router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${id}${id !== step.id ? `/node/${step.id}` : ''}`, undefined, { shallow: true })
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
                router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${id}/extend`, undefined, { shallow: true })
              }}
            />
          </Breadcrumbs>
        </div>
        <StoryProvider krg={krg} metapath={data.metapath}>
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
                    setUserPlaybook({ public: userPlaybook?.public || false })
                    setPlaybookMetadata(metadata => ({ ...metadata, id: data.metapath[data.metapath.length-1].playbook_metadata?.id || '' }))
                    setUpdateRequired(false)
                    router.push(`${session_id ? `/session/${session_id}` : ''}/report/${id}`, undefined, { shallow: true, scroll: false })
                  })
                } else {
                  deleteUserPlaybook({
                    query: { id },
                    body: {},
                  }).then(id => {
                    setUserPlaybook(undefined)
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
                    setUserPlaybook({ public: publicPlaybook })
                    setUpdateRequired(false)
                  })
                }
              }}
            />
          </Waypoint>
          {(data.metapath||[]).filter(head => head.id in cellMetadata).map(head => (
            <Waypoint
              key={`${head.id}-${cellMetadata[head.id]?.process_visible}-${cellMetadata[head.id]?.data_visible}`}
              id={head.id}
            >
              <Cell
                session_id={session_id}
                krg={krg}
                id={id}
                head={head}
                cellMetadata={cellMetadata}
                setCellMetadata={setCellMetadata}
              />
            </Waypoint>
          ))}
        </StoryProvider>
      </SessionStatus>
    </div>
  )
}