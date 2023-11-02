import React from 'react'
import type KRG from '@/core/KRG'
import dynamic from 'next/dynamic'
import { Metapath, useMetapathOutput } from '@/app/fragments/metapath'
import Head from 'next/head'
import { useStory } from '@/app/fragments/story'

const Prompt = dynamic(() => import('@/app/fragments/graph/prompt'))

export default function Cell({ session_id, krg, id, head, autoextend }: { session_id?: string, krg: KRG, id: string, head: Metapath, autoextend: boolean }) {
  const processNode = krg.getProcessNode(head.process.type)
  const { data: { output, outputNode }, error: outputError, mutate } = useMetapathOutput({ krg, head })
  const { story } = useStory()
  const [storyText, storyCitations] = React.useMemo(() => story.split('\n\n'), [story])
  const View = outputNode ? ({ output }: { output: any }) => outputNode.view(output) : undefined
  return (
    <>
      <Head>
        <title>Playbook: {processNode.meta.label}</title>
      </Head>
      <div className="flex-grow flex flex-col">
        {'prompt' in processNode ?
          <Prompt
            session_id={session_id}
            id={id}
            krg={krg}
            head={head}
            processNode={processNode}
            outputNode={outputNode}
            output={output}
            autoextend={autoextend}
          />
          : <>
          <div className="mb-4">
            <h2 className="bp5-heading">{processNode.meta.label || processNode.spec}</h2>
            <p className="prose text-justify">{storyText}</p>
            <p className="prose text-sm text-justify whitespace-pre-line">{storyCitations}</p>
          </div>
          <div className="flex-grow flex flex-col py-4">
            {outputError ? <div className="alert alert-error prose">{outputError.toString()}</div> : null}
            {!outputNode ? <div>Loading...</div>
            : <>
                {!View || output === undefined ? <div className="prose">Loading...</div>
                : output === null ? <div className="prose">Waiting for input</div>
                : <View output={output} />}
              </>}
              <button
                className="btn btn-primary"
                onClick={async (evt) => {
                  const req = await fetch(`${session_id ? `/api/socket/${session_id}` : ''}/api/db/process/${head.process.id}/output/delete`, { method: 'POST' })
                  const res = await req.text()
                  mutate()
                }}
              >Recompute</button>
          </div>
        </>}
      </div>
    </>
  )
}
