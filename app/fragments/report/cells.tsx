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
const EditableText = dynamic(() => import('@blueprintjs/core').then(({ EditableText }) => EditableText))
const Icon = dynamic(() => import('@/app/components/icon'))

export default function Cells({ krg, id }: { krg: KRG, id: string }) {
  const descriptionRef = React.useRef<HTMLTextAreaElement>(null)
  const { data: metapath, error } = useSWRImmutableSticky<Array<Metapath>>(id ? `/api/db/fpl/${id}` : undefined)
  const [metadata, setMetadata] = React.useState({ title: '', description: undefined as string | undefined, public: false })
  return (
    <div className="flex flex-col py-4 gap-2">
      <div className="flex-grow flex-shrink bp4-card p-0">
        <div className="p-3">
          <div className="flex flex-row gap-2">
            <Icon icon={start_icon} />
            <h2 className="bp4-heading">
              <EditableText
                placeholder="Playbook title"
                value={metadata.title}
                onChange={value => {setMetadata(metadata => ({ ...metadata, title: value }))}}
              />
            </h2>
          </div>
          <p
            className={`prose w-full p-1 border border-white hover:border-gray-200 cursor-text ${metadata.description !== undefined ? 'hidden' : ''}`}
            onClick={evt => {
              setMetadata(metadata => ({ ...metadata, description: evt.currentTarget?.textContent||''.replace(/\s+$/, '') }))
              setTimeout(() => {
                descriptionRef.current?.focus()
                descriptionRef.current?.click()
              }, 100)
            }}
          >
            {(metapath||[]).map(head => <Story key={head.id} krg={krg} head={head} />)}
          </p>
          <textarea
            ref={descriptionRef}
            className={`prose w-full p-1 border border-white ${metadata.description !== undefined ? '' : 'hidden'}`}
            value={metadata.description}
            onChange={evt => {
              setMetadata(metadata => ({ ...metadata, description: evt.target.value }))
            }}
          />
          {metadata.title || metadata.description !== undefined ? (
            <div className="flex flex-row gap-2 items-center">
              <button className="bp4-button bp4-intent-success">Save</button>
              <label className="bp4-control bp4-switch mb-0">
                <input
                  type="checkbox"
                />
                <span className="bp4-control-indicator"></span>
                Make workflow public
              </label>
            </div>
          ) : null}
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
