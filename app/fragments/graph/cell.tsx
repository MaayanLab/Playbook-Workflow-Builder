import React from 'react'
import type KRG from '@/core/KRG'
import { TimeoutError } from '@/spec/error'
import dynamic from 'next/dynamic'
import { Metapath, useMetapathOutput } from '@/app/fragments/metapath'
import Head from 'next/head'
import { ReactMarkdown } from 'react-markdown/lib/react-markdown'

const Prompt = dynamic(() => import('@/app/fragments/graph/prompt'))

export default function Cell({ krg, id, head, autoextend }: { krg: KRG, id: string, head: Metapath, autoextend: boolean }) {
  const processNode = krg.getProcessNode(head.process.type)
  const { data: { output, outputNode }, error: outputError, mutate } = useMetapathOutput(krg, head)
  const View = outputNode ? ({ output }: { output: any }) => outputNode.view(output) : undefined
  return (
    <>
      <Head>
        <title>Playbook: {processNode.meta.label}</title>
      </Head>
      <div className="flex-grow flex flex-col">
        {'prompt' in processNode ?
          <Prompt
            id={id}
            krg={krg}
            head={head}
            processNode={processNode}
            output={output}
            autoextend={autoextend}
          />
          : <>
          <div className="mb-4">
            <h2 className="bp4-heading">{processNode.meta.label || processNode.spec}</h2>
            {processNode.meta.description ? <p className="bp4-ui-text"><ReactMarkdown>{processNode.meta.description}</ReactMarkdown></p> : null}
          </div>
          <div className="flex-grow flex flex-col py-4">
            {outputError && !(outputError instanceof TimeoutError) ? <div className="alert alert-error">{outputError.toString()}</div> : null}
            {!outputNode ? <div>Loading...</div>
            : <>
                <h2 className="bp4-heading">{outputNode.meta.label || outputNode.spec}</h2>
                {!View || output === undefined ? <div>Loading...</div>
                : output === null ? <div>Waiting for input</div>
                : <View output={output} />}
              </>}
              <button
                className="btn btn-primary"
                onClick={async (evt) => {
                  const req = await fetch(`/api/db/process/${head.process.id}/output/delete`, { method: 'POST' })
                  const res = await req.text()
                  await mutate()
                }}
              >Recompute</button>
          </div>
        </>}
      </div>
    </>
  )
}
