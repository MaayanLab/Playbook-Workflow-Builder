import React from 'react'
import type KRG from '@/core/KRG'
import Link from 'next/link'
import { status_awaiting_input_icon, status_complete_icon, status_waiting_icon, status_alert_icon, view_in_graph_icon, fork_icon, func_icon, variable_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { Metapath, useMetapathOutput } from '@/app/fragments/metapath'
import classNames from 'classnames'
import { ReactMarkdown } from 'react-markdown/lib/react-markdown'
import { CellMetadata } from './cells'

const Prompt = dynamic(() => import('@/app/fragments/report/prompt'))
const Icon = dynamic(() => import('@/app/components/icon'))

export default function Cell({ krg, id, head, metadata, setMetadata }: { krg: KRG, id: string, head: Metapath, metadata: Record<string, CellMetadata>, setMetadata: React.Dispatch<React.SetStateAction<Record<string, CellMetadata>>> }) {
  const processNode = krg.getProcessNode(head.process.type)
  const { data: { outputNode = undefined, output = undefined } = {}, isLoading } = useMetapathOutput(krg, head)
  const View = outputNode ? outputNode.view : undefined
  return (
    <>
      {!('prompt' in processNode) ? <div className="flex-grow flex-shrink items-center overflow-auto bp4-card p-0">
        <div className="collapse collapse-arrow">
          <input type="checkbox" checked={metadata[head.id].processVisible} onChange={evt => {setMetadata((metadata) => ({ ...metadata, [head.id]: { ...metadata[head.id], processVisible: evt.target.checked } }))}} />
          <div className="collapse-title flex flex-row gap-2">
            <Icon icon={processNode.meta.icon || func_icon} />
            <h2 className="bp4-heading">
              {metadata[head.id].label ? metadata[head.id].label
                : processNode.meta.label ? processNode.meta.label
                : processNode.spec}
            </h2>
          </div>
          <div className="collapse-content">
            <p className="bp4-ui-text">
              <ReactMarkdown>
                {metadata[head.id].description ? metadata[head.id].description
                  : processNode.meta.description ? processNode.meta.description
                  : ''}
              </ReactMarkdown>
            </p>
          </div>
        </div>
        <div className={classNames('border-t-secondary border-t-2 mt-2', { 'hidden': !metadata[head.id].processVisible })}>
          <Link href={`/graph/${id}/node/${head.id}`}>
            <button className="bp4-button bp4-minimal">
              <Icon icon={view_in_graph_icon} />
            </button>
          </Link>
        </div>
      </div> : null}
      <div className="flex-grow flex-shrink items-center overflow-auto bp4-card p-0">
        {'prompt' in processNode ?
          <Prompt
            id={id}
            krg={krg}
            head={head}
            processNode={processNode}
            output={output}
          />
          : <div className="collapse collapse-arrow">
          <input type="checkbox" checked={metadata[head.id].dataVisible} onChange={evt => {setMetadata((metadata) => ({ ...metadata, [head.id]: { ...metadata[head.id], dataVisible: evt.target.checked } }))}} />
          <div className="collapse-title flex flex-row gap-2">
            <Icon icon={(outputNode && outputNode.meta.icon) || variable_icon} />
            <h2 className="bp4-heading">{(outputNode && (outputNode.meta.label || processNode.spec)) || "Loading"}</h2>
          </div>
          <div className="collapse-content">
            {outputNode && View && output ? View(output) : isLoading ? 'Waiting for results' : 'Waiting for input'}
          </div>
        </div>}
        <div className={classNames('border-t-secondary border-t-2 mt-2', { 'hidden': !metadata[head.id].dataVisible })}>
          <Link href={`/graph/${id}/node/${head.id}`}>
            <button className="bp4-button bp4-minimal">
              <Icon icon={view_in_graph_icon} />
            </button>
          </Link>
          <Link href={`/graph/${id}/node/${head.id}/extend`}>
            <button className="bp4-button bp4-minimal">
              <Icon icon={fork_icon} color="black" />
            </button>
          </Link>
          <button className="bp4-button bp4-minimal" disabled>
            {isLoading ?
              <Icon icon={status_waiting_icon} color="#DAA520" />
              : (outputNode ?
                  (output ?
                    (outputNode.spec === 'Error' ?
                      <Icon icon={status_alert_icon} color="#DC143C" />
                      : <Icon icon={status_complete_icon} color="#008000" />)
                    : <Icon icon={status_awaiting_input_icon} color="#B8860B" />)
                  : <Icon icon={status_waiting_icon} color="#DAA520" />)}
          </button>
        </div>
      </div>
    </>
  )
}
