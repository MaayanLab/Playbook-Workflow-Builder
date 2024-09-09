import React from 'react'
import type KRG from '@/core/KRG'
import Link from 'next/link'
import { status_awaiting_input_icon, status_complete_icon, status_waiting_icon, status_alert_icon, view_in_graph_icon, fork_icon, func_icon, variable_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { Metapath, useMetapathOutput } from '@/app/fragments/metapath'
import classNames from 'classnames'
import { useStory } from '@/app/fragments/story'
import { Waypoint } from '@/app/components/waypoint'
import SafeRender from '@/utils/saferender'
import { AbstractPart, FigureCaption, Methods } from './story'

const Markdown = dynamic(() => import('@/app/components/Markdown'))
const Prompt = dynamic(() => import('@/app/fragments/report/prompt'))
const Icon = dynamic(() => import('@/app/components/icon'))

export default function Cell({ session_id, krg, id, head, cellMetadata, setCellMetadata }: { session_id?: string, krg: KRG, id: string, head: Metapath, cellMetadata: Record<string, Exclude<Metapath['cell_metadata'], null>>, setCellMetadata: React.Dispatch<React.SetStateAction<Record<string, Exclude<Metapath['cell_metadata'], null>>>> }) {
  const { data: { outputNode = undefined, output = undefined } = {}, status, isLoading } = useMetapathOutput({ krg, head })
  const story = useStory()
  const processNode = krg.getProcessNode(head.process.type)
  const currentCellMetadata = cellMetadata[head.id] ?? {}
  if (!processNode) return <div className="alert alert-error">Error: {head.process.type} does not exist</div>
  return (
    <>
      <Waypoint id={`${head.id}:process`}>
        {!('prompt' in processNode) ? <div className="flex-grow flex-shrink items-center overflow-auto bp5-card p-0">
          <div className="collapse collapse-arrow text-black dark:text-white">
            <input type="checkbox" checked={currentCellMetadata.process_visible} onChange={evt => {setCellMetadata((cellMetadata) => ({ ...cellMetadata, [head.id]: { ...currentCellMetadata, process_visible: evt.target.checked, id: '' } }))}} />
            <div className="collapse-title flex flex-col gap-2">
              <div className="flex flex-row gap-2 z-10">
                <Icon icon={processNode.meta.icon || func_icon} title={processNode.meta.label} className="fill-black dark:fill-white" />
                <h2 className="bp5-heading">
                  {currentCellMetadata.label ? currentCellMetadata.label
                    : processNode.meta.label ? processNode.meta.label
                    : processNode.spec}
                </h2>
              </div>
              <AbstractPart id={head.id} story={story} />
            </div>
            <div className="collapse-content">
              <p className="bp5-ui-text">
                {currentCellMetadata.description ? <Markdown>{currentCellMetadata.description}</Markdown> : null}
                <Methods id={head.id} story={story} />
              </p>
            </div>
          </div>
          <div className={classNames('border-t-secondary border-t-2 mt-2', { 'hidden': !currentCellMetadata.process_visible })}>
            <Link href={`${session_id ? `/session/${session_id}` : ''}/graph/${id}/node/${head.id}`}>
              <button className="bp5-button bp5-minimal">
                <Icon icon={view_in_graph_icon} />
              </button>
            </Link>
          </div>
        </div> : null}
      </Waypoint>
      <Waypoint id={`${head.id}:data`}>
        <div className="flex-grow flex-shrink items-center overflow-auto bp5-card p-0">
          {'prompt' in processNode ?
            <Prompt
              session_id={session_id}
              id={id}
              krg={krg}
              head={head}
              processNode={processNode}
              outputNode={outputNode}
              output={output}
              status={status}
            />
            : <div className="collapse collapse-arrow text-black dark:text-white">
            <input type="checkbox" checked={currentCellMetadata.data_visible} onChange={evt => {setCellMetadata((cellMetadata) => ({ ...cellMetadata, [head.id]: { ...currentCellMetadata, data_visible: evt.target.checked, id: '' } }))}} />
            <div className="collapse-title flex flex-col gap-2 overflow-visible">
              <div className="flex flex-row gap-2 z-10">
                {outputNode ?
                  <>
                    <Icon icon={outputNode.meta.icon || variable_icon} title={outputNode.meta.label} className="fill-black dark:fill-white" />
                    <h2 className="bp5-heading">{outputNode.meta.label}</h2>
                  </>
                  :
                  <>
                    <Icon icon={processNode.output.meta.icon || variable_icon} title={processNode.output.meta.label} className="fill-black dark:fill-white" />
                    <h2 className="bp5-heading">{processNode.output.meta.label}</h2>
                  </>}
              </div>
            </div>
            <div className="collapse-content flex flex-col">
              {status ? (
                <div className="alert shadow-lg place-content-start align-middle">
                  <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" className="stroke-info shrink-0 w-6 h-6"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"></path></svg>
                  <code className="prose max-w-none whitespace-pre-line">{status}</code>
                </div>
              ) : null}
              {outputNode?.view && output ? <>
                <SafeRender component={outputNode.view} props={output} />
                <FigureCaption id={head.id} story={story} />
              </> : isLoading ? 'Waiting for results' : 'Waiting for input'}
            </div>
          </div>}
          <div className={classNames('border-t-secondary border-t-2 mt-2', { 'hidden': !('prompt' in processNode) && !currentCellMetadata.data_visible })}>
            <Link href={`${session_id ? `/session/${session_id}` : ''}/graph/${id}/node/${head.id}`}>
              <button className="bp5-button bp5-minimal">
                <Icon icon={view_in_graph_icon} className="fill-black dark:fill-white" />
              </button>
            </Link>
            <Link href={`${session_id ? `/session/${session_id}` : ''}/graph/${id}/node/${head.id}/extend`}>
              <button className="bp5-button bp5-minimal">
                <Icon icon={fork_icon} className="fill-black dark:fill-white" />
              </button>
            </Link>
            <button className="bp5-button bp5-minimal" disabled>
              {isLoading ?
                <Icon icon={status_waiting_icon} className="fill-yellow-500" />
                : (outputNode ?
                    (output ?
                      (outputNode.spec === 'Error' ?
                        <Icon icon={status_alert_icon} className="fill-red-500" />
                        : <Icon icon={status_complete_icon} className="fill-green-500" />)
                      : <Icon icon={status_awaiting_input_icon} className="fill-yellow-600" />)
                    : <Icon icon={status_waiting_icon} className="fill-yellow-300" />)}
            </button>
          </div>
        </div>
      </Waypoint>
    </>
  )
}
