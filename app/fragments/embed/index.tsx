import React from 'react'
import dynamic from 'next/dynamic'
import { z } from 'zod'
import { useExRouter } from '@/app/fragments/ex-router'
import { Metapath, MetapathProvider, useFPL, useKRG, useMetapathOutput } from '@/app/fragments/metapath'
import { StoryProvider } from '@/app/fragments/story'
import SafeRender from '@/utils/saferender'
import KRG from '@/core/KRG'


const ParamType = z.union([
  z.object({ session_id: z.string().optional(), graph_id: z.string(), node_id: z.string() }),
  z.object({ session_id: z.string().optional(), graph_id: z.string() }),
  z.object({ session_id: z.string().optional() }),
  z.undefined(),
])

function Cell({ krg, head }: { krg: KRG, head: Metapath }) {
  const { data: { output, outputNode }, status, error: outputError, mutate } = useMetapathOutput({ krg, head })
  if (!outputNode || !outputNode?.view || output === undefined) return <div className="prose">Loading...</div>
  if (output === null) return <div className="prose max-w-none">Waiting for input</div>
  return <><SafeRender component={outputNode.view} props={output} /><div className="embed-ready hidden"></div></>
}

function Main({ session_id, graph_id, node_id }: { session_id?: string, graph_id: string, node_id: string }) {
  const krg = useKRG()
  const { data: metapath = [] } = useFPL(graph_id)
  const head = React.useMemo(() => metapath.filter(({ id }) => id === node_id)[0], [metapath, node_id])
  return (
    <main className="flex-grow flex flex-col">
      <StoryProvider krg={krg} metapath={metapath}>
        {head ? <Cell krg={krg} head={head} /> : null}
      </StoryProvider>
    </main>
  )
}

export default function App({
  routeParams,
}: {
  routeParams?: { session_id?: string, graph_id?: string, node_id?: string },
} = {}) {
  const router = useExRouter()
  const params = ParamType.parse(routeParams ?? router.query)
  const graph_id = typeof params !== 'undefined' && 'graph_id' in params && params.graph_id ? params.graph_id : 'start'
  const node_id = typeof params !== 'undefined' && 'node_id' in params && params.node_id ? params.node_id : graph_id
  return (
    <>
      <MetapathProvider session_id={params?.session_id}>
        <Main session_id={params?.session_id} graph_id={graph_id} node_id={node_id} />
      </MetapathProvider>
    </>
  )
}
