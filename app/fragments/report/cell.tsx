import React from 'react'
import type KRG from '@/core/KRG'
import Link from 'next/link'
import { status_awaiting_input_icon, status_complete_icon, status_waiting_icon, status_alert_icon, view_in_graph_icon, fork_icon, func_icon, variable_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { Metapath, MetapathOutputs } from './metapath'

const Prompt = dynamic(() => import('@/app/fragments/report/prompt'))
const Icon = dynamic(() => import('@/app/components/icon'))

export default function Cell({ krg, id, head, metapathOutputs }: { krg: KRG, id: string, head: Metapath, metapathOutputs: MetapathOutputs }) {
  const processNode = krg.getProcessNode(head.process.type)
  const { data: { outputNode = undefined, output = undefined } = {}, isLoading } = metapathOutputs[head.process.id] || {}
  const View = outputNode ? outputNode.view : undefined
  return (
    <>
      <div key={head.id} className="flex-grow flex-shrink items-center overflow-auto bp4-card p-0">
        <div className="p-3">
          <div className="flex flex-row gap-2">
            <Icon icon={processNode.meta.icon || func_icon} />
            <h2 className="bp4-heading">{processNode.meta.label || processNode.spec}</h2>
          </div>
          {processNode.meta.description ? <p className="bp4-ui-text">{processNode.meta.description}</p> : null}
        </div>
        <div className="border-t-secondary border-t-2 mt-2">
          <Link href={`/graph/${id}/node/${head.id}`}>
            <button className="bp4-button bp4-minimal">
              <Icon icon={view_in_graph_icon} />
            </button>
          </Link>
        </div>
      </div>
      <div key={`${head.id}-output`} className="flex-grow flex-shrink items-center overflow-auto bp4-card p-0">
        {'prompt' in processNode ?
          <Prompt
            id={id}
            head={head}
            processNode={processNode}
            output={output}
          />
          : <div className="p-3">
          <div className="flex flex-row gap-2">
            <Icon icon={(outputNode && outputNode.meta.icon) || variable_icon} />
            <h2 className="bp4-heading">{(outputNode && (outputNode.meta.label || processNode.spec)) || "Loading"}</h2>
          </div>
          {outputNode && View && output ? View(output) : isLoading ? 'Waiting for results' : 'Waiting for input'}
        </div>}
        <div className="border-t-secondary border-t-2 mt-2">
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
