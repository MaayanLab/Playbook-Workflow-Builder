import React from 'react'
import type KRG from '@/core/KRG'
import dynamic from 'next/dynamic'
import { Metapath, useMetapathOutput } from '@/app/fragments/metapath'
import Head from 'next/head'
import { useStory } from '@/app/fragments/story'
import SafeRender from '@/utils/saferender'
import { Abstract, FigureCaption, Methods, References } from './story'

const Prompt = dynamic(() => import('@/app/fragments/graph/prompt'))

export default function Cell({ session_id, krg, id, head, metapath }: { session_id?: string, krg: KRG, id: string, head: Metapath, metapath: Metapath[] }) {
  const processNode = krg.getProcessNode(head.process.type)
  const { data: { output, outputNode }, status, error: outputError, mutate } = useMetapathOutput({ krg, head })
  const story = useStory()
  const astFiltered = React.useMemo(() => story.ast.filter(part => part.tags.some(tag => metapath.some(p => tag === p.id))), [story])
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
            metapath={metapath}
          />
          : <>
          <div className="mb-4">
            <h2 className="bp5-heading">{processNode.meta.label || processNode.spec}</h2>
            <Abstract story={{ ...story, ast: astFiltered }} />
            <Methods id={head.id} story={{ ...story, ast: astFiltered }} />
          </div>
          <div className="flex-grow flex flex-col py-4">
            {outputError ? <div className="alert alert-error prose max-w-none">{outputError.toString()}</div> : null}
            {status ? (
              <div className="alert shadow-lg place-content-start align-middle">
                <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" className="stroke-info shrink-0 w-6 h-6"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"></path></svg>
                <code className="prose max-w-none whitespace-pre-line">{status}</code>
              </div>
            ) : null}
            {!outputNode ? <div>Loading...</div>
              : !outputNode?.view || output === undefined ? <div className="prose">Loading...</div>
              : output === null ? <div className="prose max-w-none">Waiting for input</div>
              : <SafeRender component={outputNode.view} props={output} />}
            <FigureCaption id={head.id} story={{ ...story, ast: astFiltered }} />
            <References story={{ ...story, ast: astFiltered }} />
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
