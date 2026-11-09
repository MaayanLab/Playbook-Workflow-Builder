import { notFound } from 'next/navigation'
import Report from '@/app/app/_clients/report'
import ClientPage from '@/app/app/_clients/client-page'

export default async function SessionReportPage(props: { params: Promise<{ session_id: string, path?: string[] }> }) {
  const params = await props.params;
  if ((params.path?.length ?? 0) > 1) notFound()
  return <ClientPage><Report routeParams={{ session_id: params.session_id, id: params.path?.[0] }} /></ClientPage>
}
