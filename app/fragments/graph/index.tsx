import React from 'react'
import dynamic from 'next/dynamic'
import { z } from 'zod'
import { useExRouter } from '@/app/fragments/ex-router'
import { MetapathProvider } from '@/app/fragments/metapath'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const Graph = dynamic(() => import('@/app/fragments/graph/graph'))
const Chat = dynamic(() => import('@/app/fragments/chat/chat'))

const ParamType = z.union([
  z.object({ session_id: z.string().optional(), graph_id: z.string(), node_id: z.string(), thread_id: z.string().optional() }),
  z.object({ session_id: z.string().optional(), graph_id: z.string(), thread_id: z.string().optional() }),
  z.object({ session_id: z.string().optional(), thread_id: z.string().optional() }),
  z.undefined(),
])

export default function App({
  routeParams,
  mode,
}: {
  routeParams?: { session_id?: string, graph_id?: string, node_id?: string, thread_id?: string },
  mode?: 'extend' | 'suggest',
} = {}) {
  const router = useExRouter()
  const { params, graph_id, node_id } = React.useMemo(() => {
    const params = ParamType.parse(routeParams ?? router.query)
    const graph_id = typeof params !== 'undefined' && 'graph_id' in params && params.graph_id ? params.graph_id : 'start'
    const node_id = typeof params !== 'undefined' && 'node_id' in params && params.node_id ? params.node_id : graph_id
    return { params, graph_id, node_id }
  }, [routeParams, router.query])
  const { extend, suggest } = React.useMemo(() => ({
    extend: mode === 'extend' || router.pathname.endsWith('/extend'),
    suggest: mode === 'suggest' || router.pathname.endsWith('/suggest'),
  }), [mode, router.pathname])
  return (
    <>
      <MetapathProvider session_id={params?.session_id}>
        <Layout
          sidebar={<Chat mode="graph" graph_id={graph_id} node_id={node_id} embedded />}
        >
          <main className="grow flex flex-col mx-4">
            <Graph
              session_id={params?.session_id}
              thread_id={params?.thread_id}
              graph_id={graph_id}
              node_id={node_id}
              extend={extend}
              suggest={suggest}
            />
          </main>
        </Layout>
      </MetapathProvider>
    </>
  )
}
