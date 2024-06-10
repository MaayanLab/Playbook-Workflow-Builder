import { view_in_graph_icon } from "@/icons"
import dynamic from "next/dynamic"
import Link from "next/link"
import { useRouter } from "next/router"

const Icon = dynamic(() => import('@/app/components/icon'))

export default function GraphButton({ session_id, graph_id }: { session_id?: string, graph_id: string }) {
  const router = useRouter()
  const disabled = router.asPath.endsWith('/graph') || router.asPath.endsWith('/graph/start') || router.asPath.endsWith('/graph/extend') || router.asPath.endsWith('/graph/start/extend')
  return (
    <Link href={`${session_id ? `/session/${session_id}` : ''}/graph${graph_id === 'start' ? `/` : `/${graph_id}`}/extend`}>
      <button className='bp5-button bp5-minimal' disabled={disabled}>
        <Icon icon={view_in_graph_icon} className={disabled ? 'fill-gray-400' : 'fill-black dark:fill-white'} />
      </button>
    </Link>
  )
}
