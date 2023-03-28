import React from 'react'
import dynamic from 'next/dynamic'
import Link from 'next/link'
import { view_in_graph_icon, fork_icon, start_icon, biocompute_icon } from '@/icons'
import { useStory } from '@/app/fragments/report/story'

const ShareButton = dynamic(() => import('@/app/fragments/report/share-button'))
const BCOButton = dynamic(() => import('@/app/fragments/report/bco-button'))
const EditableText = dynamic(() => import('@blueprintjs/core').then(({ EditableText }) => EditableText))
const Icon = dynamic(() => import('@/app/components/icon'))

export default function Introduction({ id, error }: { id: string, error: any }) {
  const [metadata, setMetadata] = React.useState({ title: 'Playbook', description: undefined as string | undefined, public: false })
  const story = useStory()
  return (
    <div className="flex-grow flex-shrink bp4-card p-0">
      <div className="p-3">
        <div className="flex flex-row gap-2">
          <Icon icon={start_icon} />
          <h2 className="bp4-heading">
            <EditableText
              placeholder="Playbook title"
              value={metadata.title}
              onChange={title => {setMetadata(metadata => ({ ...metadata, title }))}}
            />
          </h2>
        </div>
        <p className="prose">
          <EditableText
            placeholder="Playbook description"
            value={metadata.description !== undefined ? metadata.description : story}
            multiline
            onChange={description => {setMetadata(metadata => ({ ...metadata, description }))}}
          />
        </p>
        {/* {metadata.title || metadata.description !== undefined ? (
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
        ) : null} */}
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
        <BCOButton id={id} metadata={metadata} />
      </div>
    </div>
  )
}