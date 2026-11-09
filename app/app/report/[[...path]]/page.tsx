import { notFound } from 'next/navigation'
import Report from '@/app/app/_clients/report'
import { workflowExists } from '@/app/app/_lib/validate-workflow'
import ClientPage from '@/app/app/_clients/client-page'

export default async function ReportPage(props: { params: Promise<{ path?: string[] }> }) {
  const params = await props.params;
  if ((params.path?.length ?? 0) > 1) notFound()
  const id = params.path?.[0]
  if (id && !(await workflowExists(id))) notFound()
  return <ClientPage><Report routeParams={{ id }} /></ClientPage>
}
