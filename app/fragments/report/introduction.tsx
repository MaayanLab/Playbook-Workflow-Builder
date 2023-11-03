import React from 'react'
import dynamic from 'next/dynamic'
import Link from 'next/link'
import Head from 'next/head'
import { view_in_graph_icon, fork_icon, start_icon, share_icon, extend_icon, func_icon, variable_icon } from '@/icons'
import { useStory } from '@/app/fragments/story'
import { useChatGPT } from '@/app/fragments/report/chatgpt'
import classNames from 'classnames'
import { Metapath } from '../metapath'
import { useRouter } from 'next/router'
import KRG from '@/core/KRG'
import * as dict from '@/utils/dict'

const SaveButton = dynamic(() => import('@/app/fragments/report/save-button'))
const LinkButton = dynamic(() => import('@/app/fragments/report/link-button'))
const BCOButton = dynamic(() => import('@/app/fragments/report/bco-button'))
const ExportButton = dynamic(() => import('@/app/fragments/report/export-button'))
const EditableText = dynamic(() => import('@blueprintjs/core').then(({ EditableText }) => EditableText))
const Icon = dynamic(() => import('@/app/components/icon'))

const Breadcrumbs = dynamic(() => import('@/app/fragments/breadcrumbs').then(({ Breadcrumbs }) => Breadcrumbs))
const DataBreadcrumb = dynamic(() => import('@/app/fragments/graph/breadcrumb').then(({ DataBreadcrumb }) => DataBreadcrumb))
const ProcessBreadcrumb = dynamic(() => import('@/app/fragments/graph/breadcrumb').then(({ ProcessBreadcrumb }) => ProcessBreadcrumb))

export default function Introduction({
  session_id,
  id,
  krg,
  metapath,
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
  const { story } = useStory()
  const [storyText, storyCitations] = React.useMemo(() => story.split('\n\n'), [story])
  const router = useRouter()
  const process_to_step = React.useMemo(() => metapath ? dict.init(metapath.map(h => ({ key: h.process.id, value: `${h.id}:${h.process.id}` }))) : {}, [metapath])
  const head = React.useMemo(() => metapath ? metapath[metapath.length - 1] : undefined, [metapath])
  return (
    <>
      <Head>
        <title>Playbook Report{playbookMetadata.title ? `: ${playbookMetadata.title}` : null}</title>
      </Head>
      <div className="sticky top-0 left-0 z-50 bg-white w-full">
        <Breadcrumbs>
          <DataBreadcrumb
            key="start"
            index={0}
            id="start"
            label="Start"
            active={false}
            icon={[start_icon]}
            parents={[]}
            onClick={() => {
              router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${id}/node/start`, undefined, { shallow: true })
            }}
          />
          {metapath.flatMap((step, i) => {
            const process = krg.getProcessNode(step.process.type)
            if (process === undefined) return []
            return [
              <ProcessBreadcrumb
                key={step.id}
                index={i * 2 + 1}
                id={step.id}
                label={process.meta.label}
                head={step}
                active={false}
                icon={process.meta.icon || [func_icon]}
                parents={dict.isEmpty(step.process.inputs) ? ['start'] : dict.values(step.process.inputs).map(({ id }) => process_to_step[id])}
                onClick={() => {
                  router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${id}${id !== step.id ? `/node/${step.id}` : ''}`, undefined, { shallow: true })
                }}
              />,
              <DataBreadcrumb
                key={`${step.id}:${step.process.id}`}
                index={i * 2 + 2}
                id={`${step.id}:${step.process.id}`}
                label={process.output.meta.label}
                head={step}
                active={false}
                icon={process.output.meta.icon || [variable_icon]}
                parents={[step.id]}
                onClick={() => {
                  router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${id}${id !== step.id ? `/node/${step.id}` : ''}`, undefined, { shallow: true })
                }}
              />,
            ]
          })}
          <ProcessBreadcrumb
            key="extend"
            index={metapath.length * 2 + 1}
            id="extend"
            label="Extend"
            active={false}
            icon={extend_icon}
            parents={[head ? `${head.id}:${head.process.id}` : `start`]}
            onClick={() => {
              router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${id}/extend`, undefined, { shallow: true })
            }}
          />
        </Breadcrumbs>
      </div>
      <div className="flex-grow flex-shrink bp5-card p-0">
        <div className="p-3">
          <div className="flex flex-row gap-2">
            <Icon icon={start_icon} className="fill-black dark:fill-white" />
            <h2 className="bp5-heading">
              <EditableText
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
                    const gpt_summary = await augmentWithChatGPT(story)
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
          <div className="prose">
            {playbookMetadata.summary === 'auto' ?
              <>
                <p className="prose-lg text-justify mt-1">{storyText}</p>
                <div className="prose-sm text-justify whitespace-pre-line">{storyCitations}</div>
              </>
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
        {error ? <div className="alert alert-error prose">{error}</div> : null}
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
                playbookMetadata.summary === 'auto' ? story
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
          />
        </div>
      </div>
    </>
  )
}