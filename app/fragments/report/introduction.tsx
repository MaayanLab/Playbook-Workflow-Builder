import React from 'react'
import dynamic from 'next/dynamic'
import Link from 'next/link'
import Head from 'next/head'
import { view_in_graph_icon, fork_icon, start_icon, share_icon, extend_icon, func_icon, variable_icon } from '@/icons'
import { useStory } from '@/app/fragments/story'
import { useChatGPT } from '@/app/fragments/report/chatgpt'
import classNames from 'classnames'
import { Metapath } from '../metapath'
import KRG from '@/core/KRG'
import { Abstract, abstractText } from './story'

const SaveButton = dynamic(() => import('@/app/fragments/report/save-button'))
const LinkButton = dynamic(() => import('@/app/fragments/report/link-button'))
const BCOButton = dynamic(() => import('@/app/fragments/report/bco-button'))
const ExportButton = dynamic(() => import('@/app/fragments/report/export-button'))
const EditableText = dynamic(() => import('@blueprintjs/core').then(({ EditableText }) => EditableText))
const Icon = dynamic(() => import('@/app/components/icon'))
const Linkify = dynamic(() => import('@/utils/linkify'))

export default function Introduction({
  session_id,
  id,
  userPlaybook,
  playbookMetadata,
  setPlaybookMetadata,
  toggleSave,
  togglePublic,
  updateRequired,
  error,
}: {
  session_id?: string,
  id: string,
  krg: KRG,
  metapath: Metapath[],
  userPlaybook?: { public: boolean },
  playbookMetadata: Exclude<Metapath['playbook_metadata'], null>, setPlaybookMetadata: React.Dispatch<React.SetStateAction<Exclude<Metapath['playbook_metadata'], null>>>,
  updateRequired: boolean,
  toggleSave: () => void,
  togglePublic: () => void,
  error: any
}) {
  const { chatGPTAvailable, augmentWithChatGPT, isAugmentingWithChatGPT, errorAugmentingWithChatGPT } = useChatGPT({ session_id })
  const story = useStory()
  return (
    <>
      <Head>
        <title>Playbook Report{playbookMetadata.title ? `: ${playbookMetadata.title}` : null}</title>
      </Head>
      <div className="flex-grow flex-shrink bp5-card p-0">
        {error ? <div className="alert alert-error prose max-w-none">{error}</div> : null}
        <div className="p-3">
          <div className="flex flex-row gap-2">
            <Icon icon={start_icon} className="fill-black dark:fill-white" />
            <h2 className="bp5-heading w-full">
              <EditableText
                multiline
                placeholder="Playbook title"
                value={playbookMetadata?.title || ''}
                onChange={title => {setPlaybookMetadata(playbookMetadata => ({ ...playbookMetadata, title, id: '' }))}}
              />
            </h2>
          </div>
          <div className="tabs">
            <button
              className={classNames('tab tab-lifted', { 'tab-active': playbookMetadata.summary === 'auto' })}
              onClick={evt => {setPlaybookMetadata(({ summary, ...playbookMetadata }) => ({ ...playbookMetadata, summary: 'auto', id: '' }))}}
            >Auto-Generated Summary</button>
            <div className="tooltip" data-tip={!chatGPTAvailable && !playbookMetadata.gpt_summary ? errorAugmentingWithChatGPT : undefined}>
              <button
                disabled={!chatGPTAvailable && !playbookMetadata.gpt_summary}
                className={classNames('tab tab-lifted', { 'tab-active': playbookMetadata.summary === 'gpt', 'cursor-not-allowed': !chatGPTAvailable && !playbookMetadata.gpt_summary })}
                onClick={async (evt) => {
                  setPlaybookMetadata(({ summary, ...playbookMetadata }) => ({ ...playbookMetadata, summary: 'gpt', id: '' }))
                  if (!playbookMetadata.gpt_summary) {
                    const gpt_summary = await augmentWithChatGPT(abstractText(story))
                    if (gpt_summary) setPlaybookMetadata(playbookMetadata => ({ ...playbookMetadata, gpt_summary, id: '' }))
                  }
                }}
              >GPT-Augmented Summary</button>
            </div>
            <button
              className={classNames('tab tab-lifted', { 'tab-active': playbookMetadata.summary === 'manual' })}
              onClick={evt => {setPlaybookMetadata(({ summary, ...playbookMetadata }) => ({ ...playbookMetadata, summary: 'manual', id: '' }))}}
            >Manual Summary</button>
          </div>
          <div className="prose max-w-none">
            {playbookMetadata.summary === 'auto' ?
              <Abstract story={story} />
            : playbookMetadata.summary === 'manual' ?
              <p className="prose-lg text-justify mt-1 whitespace-pre-line">
                <EditableText
                  placeholder="Add your manual summary here to be included when publishing."
                  value={playbookMetadata.description || ''}
                  multiline
                  onChange={description => {setPlaybookMetadata(playbookMetadata => ({ ...playbookMetadata, description, id: '' }))}}
                />
              </p>
              : playbookMetadata.summary === 'gpt' ?
              <>
                {chatGPTAvailable && isAugmentingWithChatGPT ? <progress className="progress" /> : null}
                {chatGPTAvailable && errorAugmentingWithChatGPT ? <div className="alert alert-error prose">{errorAugmentingWithChatGPT.toString()}</div> : null}
                <p className="prose-lg text-justify mt-1 whitespace-pre-line">{playbookMetadata.gpt_summary}</p>
              </>
              : null}
          </div>
        </div>
        <div className="border-t-secondary border-t-2 mt-2">
          <Link href={`${session_id ? `/session/${session_id}` : ''}/graph${id ? `/${id}/node/start` : ``}`}>
            <button className="bp5-button bp5-minimal">
              <Icon icon={view_in_graph_icon} className="fill-black dark:fill-white" />
            </button>
          </Link>
          <Link href={`${session_id ? `/session/${session_id}` : ''}/graph${id ? `/${id}/node/start/extend` : `/start/extend`}`}>
            <button className="bp5-button bp5-minimal">
              <Icon icon={fork_icon} className="fill-black dark:fill-white" />
            </button>
          </Link>
          <SaveButton
            toggleSave={toggleSave}
            userPlaybook={userPlaybook}
            updateRequired={updateRequired}
          />
          <button className="bp5-button bp5-minimal" disabled={!userPlaybook} onClick={() => {togglePublic()}}>
            <Icon
              icon={share_icon}
              className={!userPlaybook ? 'fill-gray-400' : userPlaybook.public ? 'fill-green-500' : 'fill-black dark:fill-white'}
              title={!userPlaybook ? 'Save to Share Publicly' : 'Share Publicly'}
            />
          </button>
          <BCOButton
            id={id}
            disabled={!userPlaybook}
            metadata={{
              title: playbookMetadata.title ?? 'Playbook',
              description: (
                playbookMetadata.summary === 'auto' ? abstractText(story)
                : playbookMetadata.summary === 'gpt' ? playbookMetadata.gpt_summary
                : playbookMetadata.summary === 'manual' ? playbookMetadata.description
                : undefined
              ),
            }}
          />
          <LinkButton
            id={id}
            disabled={!userPlaybook}
          />
          <ExportButton
            id={id}
            session_id={session_id}
            metadata={{
              title: playbookMetadata.title ?? 'Playbook',
              description: (
                playbookMetadata.summary === 'auto' ? abstractText(story)
                : playbookMetadata.summary === 'gpt' ? playbookMetadata.gpt_summary
                : playbookMetadata.summary === 'manual' ? playbookMetadata.description
                : undefined
              ),
            }}
          />
        </div>
      </div>
    </>
  )
}