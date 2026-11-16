import type { Metadata } from 'next'
import { notFound } from 'next/navigation'
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
export default async function SessionGraphPage(props: { params: Promise<{ session_id: string, path?: string[] }> }) {
  const route = parseGraphPath((await props.params).path?.map(decodeParam))
  if (!route) notFound()
  return null
}
