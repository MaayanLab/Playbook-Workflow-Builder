import type { Metadata } from 'next'
import { notFound } from 'next/navigation'
import Graph from '@/app/app/_clients/graph'
import { workflowNodeExists } from '@/app/app/_lib/validate-workflow'
import ClientPage from '@/app/app/_clients/client-page'
import { decodeParam } from '@/app/app/_lib/params'

function parse(encoded: string[] = []) {
  // (page params arrive url encoded, e.g. the `,` between several nodes as `%2C`)
  const path = encoded.map(decodeParam)
  let graph_id = 'start'
  let node_id: string | undefined
  let mode: 'extend' | 'suggest' | undefined
  const remaining = [...path]
  if (remaining[0] !== 'extend' && remaining[0] !== 'suggest') graph_id = remaining.shift() ?? 'start'
  if (remaining[0] === 'node') {
    remaining.shift()
    node_id = remaining.shift()
  }
  if (remaining[0] === 'extend' || remaining[0] === 'suggest') mode = remaining.shift() as typeof mode
  if (remaining.length || (path.length && !node_id && !mode && graph_id === 'start')) return undefined
  return { graph_id, node_id, mode }
}

/** Titles which don't depend on the workflow's content (steps set their own, see DocumentTitle) */
export async function generateMetadata(props: { params: Promise<{ path?: string[] }> }): Promise<Metadata> {
  const route = parse((await props.params).path)
  if (route?.mode === 'extend') return { title: 'Playbook: Extend' }
  if (route && route.graph_id === 'start' && !route.node_id && !route.mode) return { title: 'Playbook: Start' }
  return {}
}

export default async function GraphPage(props: { params: Promise<{ path?: string[] }> }) {
  const params = await props.params;
  const route = parse(params.path)
  if (!route || !(await workflowNodeExists(route.graph_id, route.node_id))) notFound()
  return <ClientPage><Graph routeParams={{ graph_id: route.graph_id, node_id: route.node_id }} mode={route.mode} /></ClientPage>
}
