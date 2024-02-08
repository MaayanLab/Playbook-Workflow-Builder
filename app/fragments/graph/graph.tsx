import React from 'react'
import { start_icon, func_icon, variable_icon, view_report_icon, Icon as IconT, extend_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { useRouter } from 'next/router'
import { useMetapath, type Metapath } from '@/app/fragments/metapath'
import useKRG from '@/app/fragments/graph/krg'
import type KRG from '@/core/KRG'
import Link from 'next/link'
import { StoryProvider } from '@/app/fragments/story'
import * as dict from '@/utils/dict'

const Breadcrumbs = dynamic(() => import('@/app/fragments/breadcrumbs'))
const Home = dynamic(() => import('@/app/fragments/playbook/home'))
const Extend = dynamic(() => import('@/app/fragments/graph/extend'))
const Suggest = dynamic(() => import('@/app/fragments/graph/suggest'))
const Cell = dynamic(() => import('@/app/fragments/graph/cell'))
const SessionStatus = dynamic(() => import('@/app/fragments/session-status'))
const ImportButton = dynamic(() => import('@/app/fragments/graph/import-button'))
const CAVATICAButton = dynamic(() => import('@/app/fragments/graph/cavatica-button'))
const RestartButton = dynamic(() => import('@/app/fragments/graph/restart-button'))
const Icon = dynamic(() => import('@/app/components/icon'))

function buildBreadcrumbGraph({
  node_id,
  metapath,
  extend,
  suggest,
  head,
  krg,
}: {
  node_id: string,
  metapath: Metapath[],
  extend: boolean,
  suggest: boolean,
  head: Metapath,
  krg: KRG,
}) {
  const graph: {
    id: string,
    kind: 'data' | 'process',
    label: string,
    color: string,
    icon: IconT,
    parents: string[],
  }[] = []
  graph.push({
    id: 'start',
    kind: 'data',
    label: 'Start',
    color: node_id === 'start' ? '#B3CFFF' : 'lightgrey',
    icon: [start_icon],
    parents: [],
  })
  // we keep this lookup to nearest process id's concrete node
  //  this is because processes with the same content have the same id regardless of position
  //  duplicate processes cause ambiguity which is trivially rectified with this strategy
  const g: Record<string, string> = {}
  for (const h of metapath) {
    const process = krg.getProcessNode(h.process.type)
    if (process !== undefined) {
      graph.push(
        {
          id: h.id,
          kind: 'process' as 'process',
          label: process.meta.label,
          color: h.id === node_id ? '#B3CFFF' : 'lightgrey',
          icon: process.meta.icon || [func_icon],
          parents: dict.isEmpty(h.process.inputs) ? ['start'] : dict.values(h.process.inputs).map(({ id }) => g[id]),
        },
        {
          id: `${h.id}:${h.process.id}`,
          kind: 'data' as 'data',
          label: process.output.meta.label,
          color: h.id === node_id ? '#B3CFFF'
            : 'prompt' in process && h.process.data?.value === undefined ? 'pink'
            : 'lightgrey',
          icon: process.output.meta.icon || [variable_icon],
          parents: [h.id],
        },
      )
      g[h.process.id] = `${h.id}:${h.process.id}`
    }
  }
  graph.push({
    id: 'extend',
    kind: 'process' as 'process',
    label: 'Extend',
    color: extend || suggest ? '#B3CFFF' : 'lightgrey',
    icon: extend_icon,
    parents: [head ? `${head.id}:${head.process.id}` : `start`],
  })
  return graph
}

function ReportButton({ session_id, graph_id }: { session_id?: string, graph_id: string }) {
  const router = useRouter()
  const disabled = router.asPath.endsWith('/graph') || router.asPath.endsWith('/graph/start') || router.asPath.endsWith('/graph/extend') || router.asPath.endsWith('/graph/start/extend')
  return (
    <Link href={`${session_id ? `/session/${session_id}` : ''}/report${graph_id === 'start' ? `/` : `/${graph_id}`}`}>
      <button className='bp5-button bp5-minimal' disabled={disabled}>
        <Icon icon={view_report_icon} className={disabled ? 'fill-gray-400' : 'fill-black dark:fill-white'} />
      </button>
    </Link>
  )
}

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
  const metapath = useMetapath().fpl ?? []
  const head = metapath.filter(({ id }) => id === node_id)[0]
  return (
    <>
      <SessionStatus session_id={session_id}>
        <div className="flex w-auto items-center justify-center">
          <Breadcrumbs
            graph={buildBreadcrumbGraph({ node_id, metapath: metapath, extend, suggest, head, krg })}
            onclick={(_evt, id) => {
              if (id === 'extend') {
                router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${graph_id}${graph_id !== node_id ? `/node/${node_id}` : ''}/extend`, undefined, { shallow: true })
              } else {
                const focus_node_id = id.split(':')[0]
                router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${graph_id}${graph_id !== focus_node_id ? `/node/${focus_node_id}` : ''}`, undefined, { shallow: true })
              }
            }}
          />
          <ImportButton session_id={session_id} />
          <CAVATICAButton session_id={session_id} />
          <RestartButton session_id={session_id} />
          <ReportButton session_id={session_id} graph_id={graph_id} />
        </div>
        <main className="flex-grow flex flex-col">
          {suggest ?
            <Suggest session_id={session_id} krg={krg} id={graph_id} head={head} />
            : extend ?
              <Extend session_id={session_id} krg={krg} id={graph_id} head={head} metapath={metapath} />
              : node_id === 'start' ?
                <Home />
                : head ?
                  <StoryProvider krg={krg} metapath={metapathToHead(metapath, head)}>
                    <Cell session_id={session_id} krg={krg} id={graph_id} head={head} autoextend />
                  </StoryProvider>
                  : null}
        </main>
      </SessionStatus>
    </>
  )
}

