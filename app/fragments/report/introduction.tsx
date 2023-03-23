import React from 'react'
import dynamic from 'next/dynamic'
import Link from 'next/link'
import { view_in_graph_icon, fork_icon, start_icon, biocompute_icon } from '@/icons'
import { useStory } from '@/app/fragments/report/story'
import { useChatGPT } from '@/app/fragments/report/chatgpt'

const ShareButton = dynamic(() => import('@/app/fragments/report/share-button'))
const EditableText = dynamic(() => import('@blueprintjs/core').then(({ EditableText }) => EditableText))
const Icon = dynamic(() => import('@/app/components/icon'))
const Bp4Popover = dynamic(() => import('@blueprintjs/popover2').then(({ Popover2 }) => Popover2))
const Bp4Menu = dynamic(() => import('@blueprintjs/core').then(({ Menu }) => Menu))
const Bp4MenuItem = dynamic(() => import('@blueprintjs/core').then(({ MenuItem }) => MenuItem))

export default function Introduction({ id, error }: { id: string, error: any }) {
  const [metadata, setMetadata] = React.useState({ title: 'Playbook', description: undefined as undefined | string, public: false })
  const { chatGPTAvailable, augmentWithChatGPT, isAugmentingWithChatGPT, errorAugmentingWithChatGPT } = useChatGPT()
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
        <p className="prose mb-2">
          <EditableText
            placeholder="Playbook description"
            value={metadata.description !== undefined ? metadata.description : story}
            multiline
            onChange={description => {setMetadata(metadata => ({ ...metadata, description }))}}
          />
          <div className="tooltip" data-tip={!chatGPTAvailable ? errorAugmentingWithChatGPT : undefined}>
            <button
              className="btn"
              disabled={!chatGPTAvailable}
              onClick={async (evt) => {
                const description = await augmentWithChatGPT(metadata.description !== undefined ? metadata.description : story)
                setMetadata(metadata => ({ ...metadata, description }))
              }}
            >Augment with ChatGPT</button>
          </div>
          {chatGPTAvailable && isAugmentingWithChatGPT ? <progress className="progress" /> : null}
          {chatGPTAvailable && errorAugmentingWithChatGPT ? <div className="alert alert-error">{errorAugmentingWithChatGPT.toString()}</div> : null}
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
        {id ? 
          <Bp4Popover
            className="cursor-pointer"
            content={
              <Bp4Menu>
                <a href={`/api/bco/${id}?metadata=${encodeURIComponent(JSON.stringify(metadata))}`} download={`${metadata.title.replace(/ /g, '-')}-bco.json`}>
                  <Bp4MenuItem
                    icon="document"
                    text="Download BCO"
                  />
                </a>
                <Bp4MenuItem
                  icon="send-to"
                  text="Draft in BioCompute Portal"
                  onClick={async (evt) => {
                    const req = await fetch(`/api/bco/${id}/draft`, { method: 'POST' })
                    const res = await req.json()
                    console.log(res)
                  }}
                />
              </Bp4Menu>
            }
            placement="bottom"
          >
            <Icon icon={biocompute_icon} color="black" title="Create BCO" />
          </Bp4Popover>
        : null}
      </div>
    </div>
  )
}