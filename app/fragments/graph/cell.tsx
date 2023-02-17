import React from 'react'
import useSWRImmutable from 'swr/immutable'
import type { Metapath } from '@/app/fragments/graph/types'
import type KRG from '@/core/KRG'
import { TimeoutError } from '@/spec/error'
import dynamic from 'next/dynamic'

const Prompt = dynamic(() => import('@/app/fragments/graph/prompt'))

export default function Cell({ krg, id, head, autoextend }: { krg: KRG, id: string, head: Metapath, autoextend: boolean }) {
  const { data: rawOutput, error: outputError } = useSWRImmutable(`/api/db/process/${head.process.id}/output`)
  const processNode = React.useMemo(() => krg.getProcessNode(head.process.type), [head])
  const outputNode = React.useMemo(() => {
    if (!processNode) return
    if (!rawOutput || outputError) return processNode.output
    return krg.getDataNode(rawOutput.type)
  }, [processNode, rawOutput, outputError])
  const { output, decodeError } = React.useMemo(() => {
    if (outputError || !outputNode) return {}
    if (!rawOutput) return { output: rawOutput }
    try {
      return { output: outputNode.codec.decode(rawOutput.value) }
    } catch (e: any) {
      console.error(e)
      return { decodeError: process.env.NODE_ENV === 'production' ? 'An unexpected error occurred.' : e.toString() }
    }
  }, [rawOutput, outputError, outputNode])
  const View = outputNode?.view
  return (
    <div className="flex-grow flex flex-col">
      {'prompt' in processNode ?
        <Prompt
          id={id}
          head={head}
          processNode={processNode}
          output={output}
          autoextend={autoextend}
        />
        : <>
        <div className="mb-4">
          <h2 className="bp4-heading">{processNode.meta.label || processNode.spec}</h2>
          {processNode.meta.description ? <p className="bp4-ui-text">{processNode.meta.description}</p> : null}
        </div>
        <div className="flex-grow flex flex-col py-4">
          {outputError && !(outputError instanceof TimeoutError) ? <div className="alert alert-error">{outputError.toString()}</div> : null}
          {!outputNode ? <div>Loading...</div>
          : <>
              <h2 className="bp4-heading">{outputNode.meta.label || outputNode.spec}</h2>
              {decodeError ? <div className="alert alert-error">{decodeError.toString()}</div>
              : !View || output === undefined ? <div>Loading...</div>
              : output === null ? <div>Waiting for input</div>
              : View(output)}
            </>}
        </div>
      </>}
    </div>
  )
}
