import { notFound } from 'next/navigation'
import { decodeParam } from '@/app/app/_lib/params'
import { parseReportPath } from '@/app/app/_lib/routes'

/** The report view is rendered by the layout (../layout), this validates the route */
export default async function SessionReportPage(props: { params: Promise<{ session_id: string, path?: string[] }> }) {
  if (!parseReportPath((await props.params).path?.map(decodeParam))) notFound()
  return null
}
