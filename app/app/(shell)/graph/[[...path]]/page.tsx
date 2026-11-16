import type { Metadata } from 'next'
import { notFound } from 'next/navigation'
import { workflowNodeExists } from '@/app/app/_lib/validate-workflow'
import { decodeParam } from '@/app/app/_lib/params'
import { parseGraphPath } from '@/app/app/_lib/routes'

/** Titles which don't depend on the workflow's content (steps set their own, see DocumentTitle) */
export async function generateMetadata(props: { params: Promise<{ path?: string[] }> }): Promise<Metadata> {
  const route = parseGraphPath((await props.params).path?.map(decodeParam))
  if (route?.mode === 'extend') return { title: 'Playbook: Extend' }
  if (route && route.graph_id === 'start' && !route.node_id && !route.mode) return { title: 'Playbook: Start' }
  return {}
}

/** The graph view is rendered by the layout (../layout), this validates the route */
export default async function GraphPage(props: { params: Promise<{ path?: string[] }> }) {
  // (page params arrive url encoded, e.g. the `,` between several nodes as `%2C`)
  const route = parseGraphPath((await props.params).path?.map(decodeParam))
  if (!route || !(await workflowNodeExists(route.graph_id, route.node_id))) notFound()
  return null
}
