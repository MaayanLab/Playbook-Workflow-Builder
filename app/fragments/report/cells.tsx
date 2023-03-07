import React from 'react'
import dynamic from 'next/dynamic'
import type KRG from '@/core/KRG'
import Link from 'next/link'
import { view_in_graph_icon, fork_icon, start_icon } from '@/icons'
import { useSWRImmutableSticky } from '@/utils/use-sticky'
import { Metapath } from './metapath'

const ShareButton = dynamic(() => import('@/app/fragments/report/share-button'))
const Cell = dynamic(() => import('@/app/fragments/report/cell'))
const Story = dynamic(() => import('@/app/fragments/report/story'))
const Icon = dynamic(() => import('@/app/components/icon'))

export default function Cells({ krg, id }: { krg: KRG, id: string }) {
  const { data: metapath, error } = useSWRImmutableSticky<Array<Metapath>>(id ? `/api/db/fpl/${id}` : undefined)
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
          <div className="prose">
            {(metapath||[]).map(head => <Story key={head.id} krg={krg} head={head} />)}
          </div>
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
      {(metapath||[]).map((head, index) => (
        <Cell key={head.id} krg={krg} id={id} head={head} defaultCollapse={index+1 !== metapath?.length} />
      ))}
    </div>
  )
}
