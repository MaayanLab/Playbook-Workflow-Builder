import React from 'react'
import dynamic from 'next/dynamic'
import Link from 'next/link'
import Head from 'next/head'
import { view_in_graph_icon, fork_icon, start_icon, biocompute_icon } from '@/icons'
import { useStory } from '@/app/fragments/report/story'
import { useChatGPT } from '@/app/fragments/report/chatgpt'
import classNames from 'classnames'
import useAsyncEffect from 'use-async-effect'

const ShareButton = dynamic(() => import('@/app/fragments/report/share-button'))
const BCOButton = dynamic(() => import('@/app/fragments/report/bco-button'))
const EditableText = dynamic(() => import('@blueprintjs/core').then(({ EditableText }) => EditableText))
const Icon = dynamic(() => import('@/app/components/icon'))

export default function Introduction({ id, error }: { id: string, error: any }) {
  const [metadata, setMetadata] = React.useState({
    title: 'Playbook',
    description: undefined as undefined | string,
    gpt_summary: undefined as undefined | string,
    summary: 'auto' as 'auto' | 'manual' | 'gpt',
    public: false,
  })
  useAsyncEffect(async (isMounted) => {
    setMetadata((metadata) => ({
      ...metadata,
      summary: metadata.summary === 'gpt' ? 'auto' : metadata.summary,
      gpt_summary: undefined,
    }))
    const published = (await import('@/app/public/playbooksDemo')).default
    if (!isMounted()) return
    published.filter((playbook) => id === playbook.id).forEach(playbook => {
      setMetadata((metadata) => ({
        ...metadata,
        title: playbook.label,
        description: playbook.description,
        gpt_summary: playbook.gpt_summary,
      }))
    })
  }, [id])
  const { chatGPTAvailable, augmentWithChatGPT, isAugmentingWithChatGPT, errorAugmentingWithChatGPT } = useChatGPT()
  const story = useStory()
  const [storyText, storyCitations] = story.split('\n\n')
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
                onChange={title => {setMetadata(metadata => ({ ...metadata, title }))}}
              />
            </h2>
          </div>
          <div className="prose">
            <div className="tabs">
              <button
                className={classNames('tab tab-lifted', { 'tab-active': metadata.summary === 'auto' })}
                onClick={evt => {setMetadata(({ summary, ...metadata }) => ({ ...metadata, summary: 'auto' }))}}
              >Auto-Generated Summary</button>
              <div className="tooltip" data-tip={!chatGPTAvailable ? errorAugmentingWithChatGPT : undefined}>
                <button
                  disabled={!chatGPTAvailable}
                  className={classNames('tab tab-lifted', { 'tab-active': metadata.summary === 'gpt', 'cursor-not-allowed': !chatGPTAvailable && !metadata.gpt_summary })}
                  onClick={async (evt) => {
                    setMetadata(({ summary, ...metadata }) => ({ ...metadata, summary: 'gpt' }))
                    if (!metadata.gpt_summary) {
                      const gpt_summary = await augmentWithChatGPT(story)
                      setMetadata((metadata) => ({ ...metadata, gpt_summary }))
                    }
                  }}
                >GPT-Augmented Summary</button>
              </div>
              <button
                className={classNames('tab tab-lifted', { 'tab-active': metadata.summary === 'manual' })}
                onClick={evt => {setMetadata(({ summary, ...metadata }) => ({ ...metadata, summary: 'manual' }))}}
              >Manual Summary</button>
            </div>
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
                  onChange={description => {setMetadata(metadata => ({ ...metadata, description }))}}
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
          <BCOButton
            id={id}
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
        </div>
      </div>
    </>
  )
}