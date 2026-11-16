import { notFound } from 'next/navigation'
import Embed from '@/app/app/_clients/embed'
import { workflowNodeExists } from '@/app/app/_lib/validate-workflow'
import ClientPage from '@/app/app/_clients/client-page'
import { decodeParam } from '@/app/app/_lib/params'

export default async function EmbedPage(props: { params: Promise<{ graph_id: string, path?: string[] }> }) {
  const params = await props.params;
  // (page params arrive url encoded)
  const graph_id = decodeParam(params.graph_id)
  const path = (params.path ?? []).map(decodeParam)
  const node_id = path.length === 2 && path[0] === 'node' ? path[1] : path.length ? undefined : graph_id
  if (!node_id || !(await workflowNodeExists(graph_id, node_id))) notFound()
  return <ClientPage><Embed routeParams={{ graph_id, node_id }} /></ClientPage>
}
