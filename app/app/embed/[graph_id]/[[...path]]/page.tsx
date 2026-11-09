import { notFound } from 'next/navigation'
import Embed from '@/app/app/_clients/embed'
import { workflowNodeExists } from '@/app/app/_lib/validate-workflow'
import ClientPage from '@/app/app/_clients/client-page'

export default async function EmbedPage(props: { params: Promise<{ graph_id: string, path?: string[] }> }) {
  const params = await props.params;
  const path = params.path ?? []
  const node_id = path.length === 2 && path[0] === 'node' ? path[1] : path.length ? undefined : params.graph_id
  if (!node_id || !(await workflowNodeExists(params.graph_id, node_id))) notFound()
  return <ClientPage><Embed routeParams={{ graph_id: params.graph_id, node_id }} /></ClientPage>
}
