import React from 'react'
import dynamic from 'next/dynamic'
import Link from 'next/link'
import Head from 'next/head'
import { view_in_graph_icon, fork_icon, start_icon, save_icon, share_icon } from '@/icons'
import { useStory } from '@/app/fragments/story'
import { useChatGPT } from '@/app/fragments/report/chatgpt'
import classNames from 'classnames'
import type { CellMetadata, Playbook } from './cells'

const LinkButton = dynamic(() => import('@/app/fragments/report/link-button'))
const BCOButton = dynamic(() => import('@/app/fragments/report/bco-button'))
const EditableText = dynamic(() => import('@blueprintjs/core').then(({ EditableText }) => EditableText))
const Icon = dynamic(() => import('@/app/components/icon'))

export default function Introduction({ id, metadata, setMetadata, playbook, setPlaybook, error }: {
  id: string,
  metadata: Record<string, CellMetadata>, setMetadata: React.Dispatch<React.SetStateAction<Record<string, CellMetadata>>>,
  playbook: Playbook, setPlaybook: React.Dispatch<React.SetStateAction<Playbook>>,
  error: any }) {
  const { chatGPTAvailable, augmentWithChatGPT, isAugmentingWithChatGPT, errorAugmentingWithChatGPT } = useChatGPT()
  const story = useStory()
  const [storyText, storyCitations] = React.useMemo(() => story.split('\n\n'), [story])
  return (
    <>
      <Head>
        <title>Playbook Report: {metadata[id].label}</title>
      </Head>
      <div className="flex-grow flex-shrink bp4-card p-0">
        <div className="p-3">
          <div className="flex flex-row gap-2">
            <Icon icon={start_icon} />
            <h2 className="bp4-heading">
              <EditableText
                placeholder="Playbook title"
                value={metadata[id].label}
                onChange={title => {setMetadata(metadata => ({ ...metadata, [id]: { ...metadata[id], title } }))}}
              />
            </h2>
          </div>
          <div className="tabs">
            <button
              className={classNames('tab tab-lifted', { 'tab-active': metadata[id].summary === 'auto' })}
              onClick={evt => {setMetadata(({ summary, ...metadata }) => ({ ...metadata, [id]: { ...metadata[id], summary: 'auto' } }))}}
            >Auto-Generated Summary</button>
            <div className="tooltip" data-tip={!chatGPTAvailable && !metadata[id].gpt_summary ? errorAugmentingWithChatGPT : undefined}>
              <button
                disabled={!chatGPTAvailable && !metadata[id].gpt_summary}
                className={classNames('tab tab-lifted', { 'tab-active': metadata[id].summary === 'gpt', 'cursor-not-allowed': !chatGPTAvailable && !metadata[id].gpt_summary })}
                onClick={async (evt) => {
                  setMetadata(({ summary, ...metadata }) => ({ ...metadata, [id]: { ...metadata[id], summary: 'gpt' } }))
                  if (!metadata[id].gpt_summary) {
                    const gpt_summary = await augmentWithChatGPT(story)
                    setMetadata((metadata) => ({ ...metadata, [id]: { ...metadata[id], gpt_summary } }))
                  }
                }}
              >GPT-Augmented Summary</button>
            </div>
            <button
              className={classNames('tab tab-lifted', { 'tab-active': metadata[id].summary === 'manual' })}
              onClick={evt => {setMetadata(({ summary, ...metadata }) => ({ ...metadata, [id]: { ...metadata[id], summary: 'manual' } }))}}
            >Manual Summary</button>
          </div>
          <div className="prose">
            {metadata[id].summary === 'auto' ?
              <>
                <p className="prose-lg mt-1">{storyText}</p>
                <div className="prose-sm whitespace-pre-line">{storyCitations}</div>
              </>
            : metadata[id].summary === 'manual' ?
              <p className="prose-lg mt-1">
                <EditableText
                  placeholder="Add your manual summary here to be included when publishing."
                  value={metadata[id].description || ''}
                  multiline
                  onChange={description => {setMetadata(metadata => ({ ...metadata, [id]: { ...metadata[id], description } }))}}
                />
              </p>
              : metadata[id].summary === 'gpt' ?
              <>
                {chatGPTAvailable && isAugmentingWithChatGPT ? <progress className="progress" /> : null}
                {chatGPTAvailable && errorAugmentingWithChatGPT ? <div className="alert alert-error">{errorAugmentingWithChatGPT.toString()}</div> : null}
                <p className="prose-lg mt-1 whitespace-pre-line">{metadata[id].gpt_summary}</p>
              </>
              : null}
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
          <button className="bp4-button bp4-minimal" onClick={() => {setPlaybook(playbook => ({ ...playbook, id: playbook.id !== id ? id : undefined, update_required: false }))}}>
            <Icon icon={save_icon} color={!playbook.id ? 'black' : playbook.update_required ? 'crimson' : 'green'} />
          </button>
          <button className="bp4-button bp4-minimal" disabled={playbook.id !== id} onClick={() => {setPlaybook(playbook => ({ ...playbook, public: playbook.public !== id ? id : undefined }))}}>
            <Icon icon={share_icon} color={playbook.id !== id ? 'gray' : playbook.public === id ? 'green' : 'black'} title="Share Publicly" />
          </button>
          <BCOButton
            id={id}
            disabled={playbook.id !== id}
            metadata={{
              title: metadata[id].label,
              description: (
                metadata[id].summary === 'auto' ? story
                : metadata[id].summary === 'gpt' ? metadata[id].gpt_summary
                : metadata[id].summary === 'manual' ? metadata[id].description
                : undefined
              ),
            }}
          />
          <LinkButton
            id={id}
            disabled={playbook.id !== id}
          />
        </div>
      </div>
    </>
  )
}