import { notFound } from 'next/navigation'
import { workflowExists } from '@/app/app/_lib/validate-workflow'
import { decodeParam } from '@/app/app/_lib/params'
import { parseReportPath } from '@/app/app/_lib/routes'

/** The report view is rendered by the layout (../layout), this validates the route */
export default async function ReportPage(props: { params: Promise<{ path?: string[] }> }) {
  // (page params arrive url encoded)
  const route = parseReportPath((await props.params).path?.map(decodeParam))
  if (!route || (route.id && !(await workflowExists(route.id)))) notFound()
  return null
}
