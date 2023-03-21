import React from 'react'
import type KRG from '@/core/KRG'
import Link from 'next/link'
import { status_awaiting_input_icon, status_complete_icon, status_waiting_icon, status_alert_icon, view_in_graph_icon, fork_icon, func_icon, variable_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { Metapath, useMetapathOutput } from '@/app/fragments/metapath'
import classNames from 'classnames'

const Prompt = dynamic(() => import('@/app/fragments/report/prompt'))
const Icon = dynamic(() => import('@/app/components/icon'))

function defaultMetadata(head: Metapath, defaultCollapse: boolean) {
  const {
    id: _,
    label,
    description,
    processVisible = false,
    dataVisible = !defaultCollapse,
  } = head.metadata ?? {}
  return { label, description, processVisible, dataVisible }
}

export default function Cell({ krg, id, head, defaultCollapse }: { krg: KRG, id: string, head: Metapath, defaultCollapse: boolean }) {
  const [metadata, setMetadata] = React.useState(defaultMetadata(head, defaultCollapse))
  const processNode = krg.getProcessNode(head.process.type)
  const { data: { outputNode = undefined, output = undefined } = {}, isLoading } = useMetapathOutput(krg, head)
  const View = outputNode ? outputNode.view : undefined
  return (
    <>
      {!('prompt' in processNode) ? <div className="flex-grow flex-shrink items-center overflow-auto bp4-card p-0">
        <div className="collapse collapse-arrow">
          <input type="checkbox" checked={metadata.processVisible} onChange={evt => {setMetadata((metadata) => ({ ...metadata, processVisible: evt.target.checked }))}} />
          <div className="collapse-title flex flex-row gap-2">
            <Icon icon={processNode.meta.icon || func_icon} />
            <h2 className="bp4-heading">
              {metadata.label ? metadata.label
                : processNode.meta.label ? processNode.meta.label
                : processNode.spec}
            </h2>
          </div>
          <div className="collapse-content">
            <p className="bp4-ui-text">
              {metadata.description ? metadata.description
                : processNode.meta.description ? processNode.meta.description
                : null}
            </p>
          </div>
        </div>
        <div className={classNames('border-t-secondary border-t-2 mt-2', { 'hidden': !metadata.processVisible })}>
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
          <input type="checkbox" checked={metadata.dataVisible} onChange={evt => {setMetadata((metadata) => ({ ...metadata, dataVisible: evt.target.checked }))}} />
          <div className="collapse-title flex flex-row gap-2">
            <Icon icon={(outputNode && outputNode.meta.icon) || variable_icon} />
            <h2 className="bp4-heading">{(outputNode && (outputNode.meta.label || processNode.spec)) || "Loading"}</h2>
          </div>
          <div className="collapse-content">
            {outputNode && View && output ? View(output) : isLoading ? 'Waiting for results' : 'Waiting for input'}
          </div>
        </div>}
        <div className={classNames('border-t-secondary border-t-2 mt-2', { 'hidden': !metadata.dataVisible })}>
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
