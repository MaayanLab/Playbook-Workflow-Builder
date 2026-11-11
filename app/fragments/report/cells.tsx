import React from 'react'
import dynamic from 'next/dynamic'
import type KRG from '@/core/KRG'
import { type Metapath, useFPL } from '@/app/fragments/metapath'
import { StoryProvider } from '@/app/fragments/story'
import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query'
import { useTRPC } from '@/app/trpc/client'
import * as dict from '@/utils/dict'
import { useExRouter } from '@/app/fragments/ex-router'
import { Breadcrumbs } from '../breadcrumbs'
import { DataBreadcrumb, ProcessBreadcrumb } from '@/app/fragments/graph/breadcrumb'
import { extend_icon, func_icon, start_icon, variable_icon } from '@/icons'
import { Waypoint, useWaypoints } from '@/app/components/waypoint'
import { References } from './story'

const Introduction = dynamic(() => import('@/app/fragments/report/introduction'))
const Cell = dynamic(() => import('@/app/fragments/report/cell'))
const SessionStatus = dynamic(() => import('@/app/fragments/session-status'))
const ShareButton = dynamic(() => import('@/app/fragments/share-button'))
const ImportButton = dynamic(() => import('@/app/fragments/graph/import-button'))
const CAVATICAButton = dynamic(() => import('@/app/fragments/graph/cavatica-button'))
const RestartButton = dynamic(() => import('@/app/fragments/graph/restart-button'))
const GraphButton = dynamic(() => import('@/app/fragments/report/graph-button'))

type ReportDraft = {
  cellMetadata: Record<string, Exclude<Metapath['cell_metadata'], null>>
  playbookMetadata: Exclude<Metapath['playbook_metadata'], null>
  cells?: Array<{
    processType: string
    metadata: Exclude<Metapath['cell_metadata'], null>
  }>
}

// App Router pages may remount when moving between graph/report URLs. Keep
// unsaved report presentation changes alive for the duration of this tab.
const reportDrafts = new Map<string, ReportDraft>()

function readReportDraft(key: string) {
  const memory = reportDrafts.get(key)
  if (memory) return memory
  if (typeof window === 'undefined') return undefined
  try {
    const value = window.sessionStorage.getItem(`pwb:report:${key}`)
    return value ? JSON.parse(value) as ReportDraft : undefined
  } catch {
    return undefined
  }
}

function writeReportDraft(key: string, draft: ReportDraft, metapath?: Metapath[]) {
  const value = metapath ? {
    ...draft,
    cells: metapath.flatMap(step => {
      const metadata = draft.cellMetadata[step.id]
      return metadata ? [{ processType: step.process.type, metadata }] : []
    }),
  } : draft
  reportDrafts.set(key, value)
  if (typeof window !== 'undefined') {
    window.sessionStorage.setItem(`pwb:report:${key}`, JSON.stringify(value))
  }
}

export default function Cells({ session_id, krg, id }: { session_id?: string, krg: KRG, id: string }) {
  const draftKey = `${session_id ?? 'local'}:${id}`
  const router = useExRouter()
  const { data: metapath } = useFPL(id)
  const trpc = useTRPC()
  const queryClient = useQueryClient()
  const { data: userPlaybook } = useQuery(trpc.playbooks.get.queryOptions({ id }))
  const mutateUserPlaybook = (value: { public: boolean } | null) => queryClient.setQueryData(trpc.playbooks.get.queryKey({ id }), value)
  const data = React.useMemo(() => metapath ? ({ metapath, userPlaybook }) : undefined, [metapath, userPlaybook])
  // in a CAVATICA session, these are handled by the worker
  const sessionContext = { trpc: { context: { session_id } } }
  const { mutateAsync: updateUserPlaybook } = useMutation(trpc.playbooks.update.mutationOptions(sessionContext))
  const { mutateAsync: publishUserPlaybook } = useMutation(trpc.playbooks.publish.mutationOptions(sessionContext))
  const { mutateAsync: deleteUserPlaybook } = useMutation(trpc.playbooks.delete.mutationOptions(sessionContext))
  const [cellMetadata, setCellMetadata] = React.useState({} as Record<string, Exclude<Metapath['cell_metadata'], null>>)
  const [playbookMetadata, setPlaybookMetadata] = React.useState({
    id: '',
    title: '',
    description: '',
    summary: 'auto',
    gpt_summary: '',
  } as Exclude<Metapath['playbook_metadata'], null>)
  const [updateRequired, setUpdateRequired] = React.useState(false)
  const previousMetapath = React.useRef<Metapath[] | undefined>(undefined)
  const pendingRebaseHead = React.useRef<string | undefined>(undefined)
  React.useEffect(() => {
    if (!data) return
    setUpdateRequired(
      playbookMetadata?.id !== data.metapath.at(-1)?.playbook_metadata?.id
      || data.metapath.some(cell => cell.cell_metadata?.id !== cellMetadata[cell.id]?.id)
    )
  }, [data, cellMetadata, playbookMetadata, userPlaybook])
  React.useLayoutEffect(() => {
    if (!metapath?.length) return
    const {
      id = '',
      title = '',
      description = '',
      summary = 'auto',
      gpt_summary = '',
    } = metapath[metapath.length-1].playbook_metadata ?? {}
    const draft = readReportDraft(draftKey)
    setPlaybookMetadata(draft?.playbookMetadata ?? { id, title, description, summary, gpt_summary })
    // Snapshot this before enqueueing the state update. Reading the ref from
    // inside the updater races with the assignment below and can make the new
    // metapath masquerade as the previous one.
    const oldMetapath = previousMetapath.current
    const oldStepIds = new Set(oldMetapath?.map(step => step.id))
    const oldProcessIds = new Set(oldMetapath?.map(step => step.process.id))
    const isExpectedRebase = pendingRebaseHead.current === metapath[metapath.length - 1]?.id
    const isRebase = oldMetapath !== undefined && (isExpectedRebase || metapath.some(step =>
      oldStepIds.has(step.id) || oldProcessIds.has(step.process.id)
    ))
    const oldStepsBySlot = new Map<string, Metapath>()
    const oldTypeCounts = new Map<string, number>()
    oldMetapath?.forEach(step => {
      const occurrence = oldTypeCounts.get(step.process.type) ?? 0
      oldTypeCounts.set(step.process.type, occurrence + 1)
      oldStepsBySlot.set(`${step.process.type}:${occurrence}`, step)
    })
    setCellMetadata(previous => {
      // A rebase can add or remove downstream steps, so equality of the whole
      // path is too strict. An unchanged step is enough to establish that this
      // is another revision of the graph currently on screen.
      const typeCounts = new Map<string, number>()
      const next = dict.init(metapath.map((element, index) => {
        const {
          id = '',
          label = '',
          description = '',
          process_visible = false,
          data_visible = index+1 !== metapath.length-1,
        } = element.cell_metadata ?? {}
        const serverValue = { id, label, description, process_visible, data_visible }
        const occurrence = typeCounts.get(element.process.type) ?? 0
        typeCounts.set(element.process.type, occurrence + 1)
        const oldStep = isRebase ? oldStepsBySlot.get(`${element.process.type}:${occurrence}`) : undefined
        const localValue = oldStep ? previous[oldStep.id] : undefined
        // An empty metadata id marks unsaved presentation changes. Carry them
        // across a content-addressed input rebase, whose FPL ids all change.
        return { key: element.id, value: localValue?.id === '' ? localValue : serverValue }
      }))
      if (!draft) return next
      const draftTypeCounts = new Map<string, number>()
      const draftBySlot = new Map(draft.cells?.map(cell => {
        const occurrence = draftTypeCounts.get(cell.processType) ?? 0
        draftTypeCounts.set(cell.processType, occurrence + 1)
        return [`${cell.processType}:${occurrence}`, cell.metadata] as const
      }))
      const currentTypeCounts = new Map<string, number>()
      return dict.init(metapath.map(step => {
        const occurrence = currentTypeCounts.get(step.process.type) ?? 0
        currentTypeCounts.set(step.process.type, occurrence + 1)
        return {
          key: step.id,
          value: draft.cellMetadata[step.id]
            ?? draftBySlot.get(`${step.process.type}:${occurrence}`)
            ?? next[step.id],
        }
      }))
    })
    if (pendingRebaseHead.current === metapath[metapath.length - 1]?.id) {
      pendingRebaseHead.current = undefined
    }
    previousMetapath.current = metapath
  }, [metapath, draftKey])
  React.useLayoutEffect(() => {
    // Do not overwrite the incoming draft with the previous render's metadata
    // while the initialization effect is remapping the new workflow IDs.
    if (!metapath?.length || metapath[metapath.length - 1].id !== id
      || metapath.some(step => !(step.id in cellMetadata))) return
    writeReportDraft(draftKey, { cellMetadata, playbookMetadata }, metapath)
  }, [draftKey, metapath, cellMetadata, playbookMetadata])
  const process_to_step = React.useMemo(() => metapath ? dict.init(metapath.map(h => ({ key: h.process.id, value: `${h.id}:${h.process.id}` }))) : {}, [metapath])
  const head = React.useMemo(() => metapath ? metapath[metapath.length - 1] : undefined, [metapath])
  const { waypoints, scrollTo } = useWaypoints()
  if (!data || !playbookMetadata || !metapath) return null
  return (
    <div className="flex flex-col py-4 gap-2">
      <SessionStatus session_id={session_id}>
        <StoryProvider krg={krg} metapath={data.metapath}>
          <Waypoint id="head" className="sticky top-0 z-20 bg-white dark:bg-current w-full flex flex-row place-items-center">
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
                  router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${id}/extend`, undefined, { shallow: true })
                }}
              />
            </Breadcrumbs>
            <ShareButton disabled={!!session_id} />
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
              userPlaybook={userPlaybook ?? undefined}
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
                    id,
                    user_playbook: { public: userPlaybook?.public || false },
                    playbook_metadata: playbookMetadataId ? { id: playbookMetadataId } : playbook_metadata,
                    cell_metadata,
                  }).then(id => {
                    mutateUserPlaybook({ public: userPlaybook?.public || false })
                    setPlaybookMetadata(metadata => ({ ...metadata, id: data.metapath[data.metapath.length-1].playbook_metadata?.id || '' }))
                    setUpdateRequired(false)
                    router.push(`${session_id ? `/session/${session_id}` : ''}/report/${id}`, undefined, { shallow: true, scroll: false })
                  })
                } else {
                  deleteUserPlaybook({ id }).then(() => {
                    mutateUserPlaybook(null)
                    setUpdateRequired(false)
                  })
                }
              }}
              togglePublic={() => {
                if (!updateRequired && userPlaybook) {
                  const publicPlaybook = !userPlaybook.public
                  publishUserPlaybook({ id, public: publicPlaybook }).then(() => {
                    mutateUserPlaybook({ public: publicPlaybook })
                    setUpdateRequired(false)
                  })
                }
              }}
            />
          </Waypoint>
          {(data.metapath||[]).map((head, index) => (
            <Cell
              key={`${head.process.type}:${index}`}
              session_id={session_id}
              krg={krg}
              id={id}
              head={head}
              collapseKey={`${session_id ?? 'local'}:${head.process.type}:${index}`}
              cellMetadata={cellMetadata}
              setCellMetadata={setCellMetadata}
              onRebase={headId => {
                pendingRebaseHead.current = headId
                writeReportDraft(`${session_id ?? 'local'}:${headId}`, { cellMetadata, playbookMetadata }, metapath)
              }}
            />
          ))}
          <References />
        </StoryProvider>
      </SessionStatus>
    </div>
  )
}
