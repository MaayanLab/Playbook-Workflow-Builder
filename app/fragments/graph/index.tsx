import React from 'react'
import dynamic from 'next/dynamic'
import { usePathname, useSearchParams } from 'next/navigation'
import { parsePathname } from '@/app/app/_lib/routes'

const Graph = dynamic(() => import('@/app/fragments/graph/graph'))

/**
 * The graph view of the route (rendered by the graph layouts, so that it persists, and keeps its
 *  state, while navigating between graphs & nodes); the shell & playbook state are provided by
 *  @/app/app/(shell)/layout.
 */
export default function App() {
  const pathname = usePathname() ?? ''
  const searchParams = useSearchParams()
  const route = React.useMemo(() => parsePathname(pathname), [pathname])
  if (route.kind !== 'graph') return null
  return (
    <main className="grow flex flex-col mx-4">
      <Graph
        session_id={route.session_id}
        thread_id={searchParams?.get('thread_id') ?? undefined}
        graph_id={route.graph_id}
        node_id={route.node_id ?? route.graph_id}
        extend={route.mode === 'extend'}
        suggest={route.mode === 'suggest'}
      />
    </main>
  )
}
