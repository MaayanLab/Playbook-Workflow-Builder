import React from 'react'
import dynamic from 'next/dynamic'
import type { FPL } from '@/core/FPPRG'
import type KRG from '@/core/KRG'
import Link from 'next/link'
import FPL2Text from '@/core/fpl2text'
import { view_in_graph_icon, fork_icon, start_icon } from '@/icons'
import { useSWRImmutableSticky } from '@/utils/use-sticky'

const ShareButton = dynamic(() => import('@/app/fragments/report/share-button'))
const Cell = dynamic(() => import('@/app/fragments/report/cell'))
const Icon = dynamic(() => import('@/app/components/icon'))

type Metapath = ReturnType<FPL['toJSON']>

export default function Cells({ krg, id }: { krg: KRG, id: string }) {
  const { data: raw_metapath, error } = useSWRImmutableSticky<Array<Metapath>>(id ? `/api/db/fpl/${id}` : undefined)
  const metapath = raw_metapath||[]
  const metapath_filtered = React.useMemo(() => {
    return metapath
      .filter((head, index) =>
        'prompt' in krg.getProcessNode(head.process.type) // any prompt nodes
        || index === metapath.length-1 // last step
      )
  }, [metapath])
  const description = React.useMemo(() => FPL2Text(krg, metapath), [krg, metapath])
  return (
    <div className="flex flex-col py-4 gap-2">
      <div className="flex-grow flex-shrink bp4-card p-0">
        <div className="p-3">
          <div className="flex flex-row gap-2">
            <Icon icon={start_icon} />
            <h2 className="bp4-heading">
              Playbook
            </h2>
          </div>
          <h2 className="prose whitespace-pre-line">
            {description}
          </h2>
        </div>
        {error ? <div className="alert alert-error">{error}</div> : null}
        <div className="border-t-secondary border-t-2 mt-2">
          <Link href={`/graph${id ? `/${id}/node/start` : ``}`}>
            <button className="bp4-button bp4-minimal">
              <Icon icon={view_in_graph_icon} />
            </button>
          </Link>
          <Link href={`/graph${id ? `/${id}/node/start/extend` : `/start/extend`}`}>
            <button className="bp4-button bp4-minimal">
              <Icon icon={fork_icon} color="black" />
            </button>
          </Link>
          <ShareButton id={id} />
        </div>
      </div>
      {metapath_filtered.map((head, index) =>
        <Cell key={index} krg={krg} index={index} id={id} head={head} />
      )}
    </div>
  )
}
