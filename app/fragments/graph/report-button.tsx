import { useRouter } from 'next/router'
import dynamic from 'next/dynamic'
import Link from 'next/link'
import { view_report_icon } from '@/icons'
const Icon = dynamic(() => import('@/app/components/icon'))

export default function ReportButton({ session_id, graph_id, thread_id }: { session_id?: string, graph_id: string, thread_id?: string }) {
  const router = useRouter()
  const disabled = router.asPath.endsWith('/graph') || router.asPath.endsWith('/graph/start') || router.asPath.endsWith('/graph/extend') || router.asPath.endsWith('/graph/start/extend') || graph_id === 'start'
  return (
    <Link href={`${session_id ? `/session/${session_id}` : ''}/report${graph_id === 'start' ? `/` : `/${graph_id}`}${thread_id ? `?thread=${thread_id}` : ''}`}>
      <button className='bp5-button bp5-minimal' disabled={disabled}>
        <Icon icon={view_report_icon} className={disabled ? 'fill-gray-400' : 'fill-black dark:fill-white'} />
      </button>
    </Link>
  )
}
