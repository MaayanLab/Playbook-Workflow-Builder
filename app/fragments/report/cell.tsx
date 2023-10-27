import React from 'react'
import type KRG from '@/core/KRG'
import Link from 'next/link'
import { status_awaiting_input_icon, status_complete_icon, status_waiting_icon, status_alert_icon, view_in_graph_icon, fork_icon, func_icon, variable_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { Metapath, useMetapathOutput } from '@/app/fragments/metapath'
import classNames from 'classnames'
import { ReactMarkdown } from 'react-markdown/lib/react-markdown'

const Prompt = dynamic(() => import('@/app/fragments/report/prompt'))
const Icon = dynamic(() => import('@/app/components/icon'))

export default function Cell({ session_id, krg, id, head, cellMetadata, setCellMetadata }: { session_id?: string, krg: KRG, id: string, head: Metapath, cellMetadata: Record<string, Exclude<Metapath['cell_metadata'], null>>, setCellMetadata: React.Dispatch<React.SetStateAction<Record<string, Exclude<Metapath['cell_metadata'], null>>>> }) {
  const { data: { outputNode = undefined, output = undefined } = {}, isLoading } = useMetapathOutput({ session_id, krg, head })
  const View = outputNode ? outputNode.view : undefined
  const processNode = krg.getProcessNode(head.process.type)
  if (!processNode) return <div className="alert alert-error">Error: {head.process.type} does not exist</div>
  return (
    <>
      {!('prompt' in processNode) ? <div className="flex-grow flex-shrink items-center overflow-auto bp5-card p-0">
        <div className="collapse collapse-arrow text-black dark:text-white">
          <input type="checkbox" checked={cellMetadata[head.id].process_visible} onChange={evt => {setCellMetadata((cellMetadata) => ({ ...cellMetadata, [head.id]: { ...cellMetadata[head.id], process_visible: evt.target.checked, id: '' } }))}} />
          <div className="collapse-title flex flex-row gap-2">
            <Icon icon={processNode.meta.icon || func_icon} className="fill-black dark:fill-white" />
            <h2 className="bp5-heading">
              {cellMetadata[head.id].label ? cellMetadata[head.id].label
                : processNode.meta.label ? processNode.meta.label
                : processNode.spec}
            </h2>
          </div>
          <div className="collapse-content">
            <p className="bp5-ui-text">
              <ReactMarkdown>
                {cellMetadata[head.id].description ?? processNode.meta.description ?? ''}
              </ReactMarkdown>
            </p>
          </div>
        </div>
        <div className={classNames('border-t-secondary border-t-2 mt-2', { 'hidden': !cellMetadata[head.id].process_visible })}>
          <Link href={`${session_id ? `/session/${session_id}` : ''}/graph/${id}/node/${head.id}`}>
            <button className="bp5-button bp5-minimal">
              <Icon icon={view_in_graph_icon} />
            </button>
          </Link>
        </div>
      </div> : null}
      <div className="flex-grow flex-shrink items-center overflow-auto bp5-card p-0">
        {'prompt' in processNode ?
          <Prompt
            session_id={session_id}
            id={id}
            krg={krg}
            head={head}
            processNode={processNode}
            output={output}
          />
          : <div className="collapse collapse-arrow text-black dark:text-white">
          <input type="checkbox" checked={cellMetadata[head.id].data_visible} onChange={evt => {setCellMetadata((cellMetadata) => ({ ...cellMetadata, [head.id]: { ...cellMetadata[head.id], data_visible: evt.target.checked, id: '' } }))}} />
          <div className="collapse-title flex flex-row gap-2">
            <Icon icon={(outputNode && outputNode.meta.icon) || variable_icon} className="fill-black dark:fill-white" />
            <h2 className="bp5-heading">{(outputNode && (outputNode.meta.label || processNode.spec)) || "Loading"}</h2>
          </div>
          <div className="collapse-content">
            {outputNode && View && output ? View(output) : isLoading ? 'Waiting for results' : 'Waiting for input'}
          </div>
        </div>}
        <div className={classNames('border-t-secondary border-t-2 mt-2', { 'hidden': !('prompt' in processNode) && !cellMetadata[head.id].data_visible })}>
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
    </>
  )
}
