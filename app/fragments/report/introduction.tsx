import React from 'react'
import dynamic from 'next/dynamic'
import Link from 'next/link'
import Head from 'next/head'
import { view_in_graph_icon, fork_icon, start_icon, save_icon, share_icon } from '@/icons'
import { useStory } from '@/app/fragments/report/story'
import { useChatGPT } from '@/app/fragments/report/chatgpt'
import classNames from 'classnames'
import type { ReportMetadata } from './cells'

const LinkButton = dynamic(() => import('@/app/fragments/report/link-button'))
const BCOButton = dynamic(() => import('@/app/fragments/report/bco-button'))
const EditableText = dynamic(() => import('@blueprintjs/core').then(({ EditableText }) => EditableText))
const Icon = dynamic(() => import('@/app/components/icon'))

export default function Introduction({ id, metadata, setMetadata, error }: { id: string, metadata: ReportMetadata, setMetadata: React.Dispatch<React.SetStateAction<ReportMetadata>>, error: any }) {
  const { chatGPTAvailable, augmentWithChatGPT, isAugmentingWithChatGPT, errorAugmentingWithChatGPT } = useChatGPT()
  const story = useStory()
  const [storyText, storyCitations] = React.useMemo(() => story.split('\n\n'), [story])
  return (
    <>
      <Head>
        <title>Playbook Report: {metadata.title}</title>
      </Head>
      <div className="flex-grow flex-shrink bp4-card p-0">
        <div className="p-3">
          <div className="flex flex-row gap-2">
            <Icon icon={start_icon} />
            <h2 className="bp4-heading">
              <EditableText
                placeholder="Playbook title"
                value={metadata.title}
                onChange={title => {setMetadata(metadata => ({ ...metadata, title, saved: metadata.saved === 'yes' ? 'pending' : 'no' }))}}
              />
            </h2>
          </div>
          <div className="tabs">
            <button
              className={classNames('tab tab-lifted', { 'tab-active': metadata.summary === 'auto' })}
              onClick={evt => {setMetadata(({ summary, ...metadata }) => ({ ...metadata, summary: 'auto', saved: metadata.saved === 'yes' ? 'pending' : 'no' }))}}
            >Auto-Generated Summary</button>
            <div className="tooltip" data-tip={!chatGPTAvailable && !metadata.gpt_summary ? errorAugmentingWithChatGPT : undefined}>
              <button
                disabled={!chatGPTAvailable && !metadata.gpt_summary}
                className={classNames('tab tab-lifted', { 'tab-active': metadata.summary === 'gpt', 'cursor-not-allowed': !chatGPTAvailable && !metadata.gpt_summary })}
                onClick={async (evt) => {
                  setMetadata(({ summary, ...metadata }) => ({ ...metadata, summary: 'gpt', saved: metadata.saved === 'yes' ? 'pending' : 'no' }))
                  if (!metadata.gpt_summary) {
                    const gpt_summary = await augmentWithChatGPT(story)
                    setMetadata((metadata) => ({ ...metadata, gpt_summary, saved: metadata.saved === 'yes' ? 'pending' : 'no' }))
                  }
                }}
              >GPT-Augmented Summary</button>
            </div>
            <button
              className={classNames('tab tab-lifted', { 'tab-active': metadata.summary === 'manual' })}
              onClick={evt => {setMetadata(({ summary, ...metadata }) => ({ ...metadata, summary: 'manual', saved: metadata.saved === 'yes' ? 'pending' : 'no' }))}}
            >Manual Summary</button>
          </div>
          <div className="prose">
            {metadata.summary === 'auto' ?
              <>
                <p className="prose-lg mt-1">{storyText}</p>
                <div className="prose-sm whitespace-pre-line">{storyCitations}</div>
              </>
            : metadata.summary === 'manual' ?
              <p className="prose-lg mt-1">
                <EditableText
                  placeholder="Add your manual summary here to be included when publishing."
                  value={metadata.description || ''}
                  multiline
                  onChange={description => {setMetadata(metadata => ({ ...metadata, description, saved: metadata.saved === 'yes' ? 'pending' : 'no' }))}}
                />
              </p>
              : metadata.summary === 'gpt' ?
              <>
                {chatGPTAvailable && isAugmentingWithChatGPT ? <progress className="progress" /> : null}
                {chatGPTAvailable && errorAugmentingWithChatGPT ? <div className="alert alert-error">{errorAugmentingWithChatGPT.toString()}</div> : null}
                <p className="prose-lg mt-1 whitespace-pre-line">{metadata.gpt_summary}</p>
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
          <button className="bp4-button bp4-minimal" onClick={() => {setMetadata(metadata => ({ ...metadata, saved: metadata.saved === 'yes' ? 'no' : 'yes' }))}}>
            <Icon icon={save_icon} color={metadata.saved === 'yes' ? 'green' : metadata.saved === 'pending' ? 'orange': 'black'} />
          </button>
          <button className="bp4-button bp4-minimal" disabled={metadata.saved !== 'yes'} onClick={() => {setMetadata(metadata => ({ ...metadata, public: !metadata.public }))}}>
            <Icon icon={share_icon} color={metadata.saved === 'yes' ? metadata.public ? 'green' : 'black' : 'gray'} title="Share Publicly" />
          </button>
          <BCOButton
            id={id}
            disabled={metadata.saved !== 'yes'}
            metadata={{
              title: metadata.title,
              description: (
                metadata.summary === 'auto' ? story
                : metadata.summary === 'gpt' ? metadata.gpt_summary
                : metadata.summary === 'manual' ? metadata.description
                : undefined
              ),
            }}
          />
          <LinkButton
            id={id}
            disabled={metadata.saved !== 'yes'}
          />
        </div>
      </div>
    </>
  )
}