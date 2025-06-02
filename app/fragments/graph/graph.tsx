import React from 'react'
import { start_icon, func_icon, variable_icon, extend_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { useRouter } from 'next/router'
import { type Metapath, useFPL } from '@/app/fragments/metapath'
import useKRG from '@/app/fragments/graph/krg'
import { StoryProvider } from '@/app/fragments/story'
import * as dict from '@/utils/dict'
import { Set } from 'immutable'

const Breadcrumbs = dynamic(() => import('@/app/fragments/breadcrumbs').then(({ Breadcrumbs }) => Breadcrumbs))
const DataBreadcrumb = dynamic(() => import('@/app/fragments/graph/breadcrumb').then(({ DataBreadcrumb }) => DataBreadcrumb))
const ProcessBreadcrumb = dynamic(() => import('@/app/fragments/graph/breadcrumb').then(({ ProcessBreadcrumb }) => ProcessBreadcrumb))
const Home = dynamic(() => import('@/app/fragments/playbook/home'))
const Extend = dynamic(() => import('@/app/fragments/graph/extend'))
const Suggest = dynamic(() => import('@/app/fragments/graph/suggest'))
const Cell = dynamic(() => import('@/app/fragments/graph/cell'))
const SessionStatus = dynamic(() => import('@/app/fragments/session-status'))
const ImportButton = dynamic(() => import('@/app/fragments/graph/import-button'))
const CAVATICAButton = dynamic(() => import('@/app/fragments/graph/cavatica-button'))
const RestartButton = dynamic(() => import('@/app/fragments/graph/restart-button'))
const ReportButton = dynamic(() => import('@/app/fragments/graph/report-button'))
const ShareButton = dynamic(() => import('@/app/fragments/share-button'))


/**
 * Find the metapath to the current head, excluding irrelevant steps
 */
function metapathToHead(metapath: Array<Metapath>, head: Metapath) {
  const relevant: Record<string, boolean> = {}
  const processMetapathMapping = dict.init(metapath.map(h => ({ key: h.process.id, value: h })))
  const Q = [head]
  while (true) {
    const el = Q.pop()
    if (!el) break
    relevant[el.id] = true
    Object.values(el.process.inputs).map(dep => {
      Q.push(processMetapathMapping[dep.id])
    })
  }
  return metapath.filter(({ id }) => relevant[id] === true)
}

export default function Graph({ session_id, graph_id, node_id, extend, suggest }: { session_id?: string, graph_id: string, node_id: string, extend: boolean, suggest: boolean }) {
  const router = useRouter()
  const krg = useKRG({ session_id })
  const { data: metapath = [] } = useFPL(graph_id)
  const node_ids = React.useMemo(() => Set(node_id.split(',')), [node_id])
  const heads = React.useMemo(() => {
    return metapath.filter(({ id }) => node_ids.has(id))
  }, [metapath, node_ids])
  const process_to_step = React.useMemo(() => dict.init(metapath.map(h => ({ key: h.process.id, value: `${h.id}:${h.process.id}` }))), [metapath])
  return (
    <>
      <SessionStatus session_id={session_id}>
        <div className="flex w-auto items-center justify-center">
          <Breadcrumbs>
            <DataBreadcrumb
              key="start"
              index={0}
              id="start"
              label="Start"
              active={node_id === 'start' && !extend}
              icon={[start_icon]}
              parents={[]}
              onClick={() => {
                router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${graph_id}${graph_id !== 'start' ? `/node/start` : ''}`, undefined, { shallow: true })
              }}
            />
            {metapath.flatMap((step, i) => {
              const process = krg.getProcessNode(step.process.type)
              if (process === undefined) return []
              return [
                <ProcessBreadcrumb
                  key={step.id}
                  index={i*2+1}
                  id={step.id}
                  label={process.meta.label}
                  head={step}
                  active={false}
                  icon={process.meta.icon || [func_icon]}
                  parents={dict.isEmpty(step.process.inputs) ? ['start'] : dict.values(step.process.inputs).map(({ id }) => process_to_step[id])}
                  onClick={(evt) => {
                    if (evt.shiftKey) {
                      const new_node_id = node_ids.has(step.id) ? node_ids.delete(step.id).join(',') : node_ids.add(step.id).join(',')
                      if (new_node_id) {
                        router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${graph_id}/node/${new_node_id}/extend`, undefined, { shallow: true })
                      }
                    } else {
                      router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${graph_id}${graph_id !== step.id ? `/node/${step.id}` : ''}`, undefined, { shallow: true })
                    }
                  }}
                />,
                <DataBreadcrumb
                  key={`${step.id}:${step.process.id}`}
                  index={i*2+2}
                  id={`${step.id}:${step.process.id}`}
                  label={process.output.meta.label}
                  head={step}
                  active={node_ids.has(step.id)}
                  icon={process.output.meta.icon || [variable_icon]}
                  parents={[step.id]}
                  onClick={(evt) => {
                    if (evt.shiftKey) {
                      const new_node_id = node_ids.has(step.id) ? node_ids.delete(step.id).join(',') : node_ids.add(step.id).join(',')
                      if (new_node_id) {
                        router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${graph_id}/node/${new_node_id}/extend`, undefined, { shallow: true })
                      }
                    } else {
                      router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${graph_id}${graph_id !== step.id ? `/node/${step.id}` : ''}`, undefined, { shallow: true })
                    }
                  }}
                />,
              ]
            })}
            <ProcessBreadcrumb
              key="extend"
              index={metapath.length*2+1}
              id="extend"
              label="Extend"
              active={extend || suggest}
              icon={extend_icon}
              parents={[
                ...(heads.length === 0 ? ['start'] : []),
                ...(heads.length > 0 ? heads.map(head => `${head.id}:${head.process.id}`) : []),
              ]}
              onClick={() => {
                router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${graph_id}${graph_id !== node_id ? `/node/${node_id}` : ''}/extend`, undefined, { shallow: true })
              }}
            />
          </Breadcrumbs>
          <ShareButton disabled={!!session_id} />
          <ImportButton session_id={session_id} />
          <CAVATICAButton session_id={session_id} />
          <RestartButton session_id={session_id} />
          <ReportButton session_id={session_id} graph_id={graph_id} />
        </div>
        <main className="flex-grow flex flex-col">
          {suggest ?
            <Suggest session_id={session_id} krg={krg} id={graph_id} head={heads[0]} />
            : extend ?
              <Extend session_id={session_id} krg={krg} id={graph_id} heads={heads} metapath={metapath} />
              : node_id === 'start' ?
                <Home />
                : heads.length === 1 ?
                  <StoryProvider krg={krg} metapath={metapath}>
                    <Cell session_id={session_id} krg={krg} id={graph_id} head={heads[0]} metapath={metapathToHead(metapath, heads[0])} />
                  </StoryProvider>
                  : null}
        </main>
      </SessionStatus>
    </>
  )
}

