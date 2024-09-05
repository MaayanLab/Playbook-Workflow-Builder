import React from 'react'
import { z } from 'zod'
import { useRouter } from 'next/router'
import { DataMetaNode, PromptMetaNode } from '@/spec/metanode'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { func_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { Metapath, useMetapathInputs } from '@/app/fragments/metapath'
import type KRG from '@/core/KRG'
import { useStory } from '@/app/fragments/story'
import classNames from 'classnames'
import { AbstractPart, FigureCaption } from './story'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function Prompt({ session_id, krg, processNode, outputNode, output, id, head, status }: { session_id?: string, krg: KRG, processNode: PromptMetaNode, outputNode?: DataMetaNode, output: any, id: string, head: Metapath, status: string | undefined }) {
  const router = useRouter()
  const { data: inputs, error } = useMetapathInputs({ krg, head })
  const story = useStory()
  const data = React.useMemo(() => {
    // invalidate data if it no longer matches the codec
    //  which could happen if an older verison of the metanode
    //  was used originally
    try {
      if (head.process.data !== null) {
        return processNode.codec.decode(head.process.data.value)
      }
    } catch (e) {}
  }, [head])
  const Component = processNode.prompt
  return (
    <div className={classNames('collapse collapse-arrow text-black dark:text-white', { 'bg-yellow-300 dark:bg-yellow-600': output === null, 'bg-green-400 dark:bg-green-700': output !== null })}>
      <input type="checkbox" defaultChecked={true} />
      <div className="collapse-title flex flex-col gap-2">
        <div className="flex flex-row gap-2 z-10">
          <Icon icon={processNode.meta.icon || func_icon} title={processNode.meta.label} className="fill-black dark:fill-white" />
          <h2 className="bp5-heading">{processNode.meta.label || processNode.spec}</h2>
        </div>
        <AbstractPart id={head.id} story={story} />
      </div>
      {error ? <div className="alert alert-error prose max-w-none">{error.toString()}</div> : null}
      <div className="collapse-content">
        {inputs !== undefined && array.intersection(dict.keys(processNode.inputs), dict.keys(inputs)).length === dict.keys(processNode.inputs).length ?
          <Component
            session_id={session_id}
            data={data}
            inputs={inputs}
            output={outputNode?.spec === processNode.output.spec ? output : undefined}
            submit={async (data) => {
              const req = await fetch(`${session_id ? `/api/socket/${session_id}` : ''}/api/db/fpl/${id}/rebase/${head.process.id}`, {
                headers: {
                  'Content-Type': 'application/json',
                },
                method: 'POST',
                body: JSON.stringify({
                  type: head.process.type,
                  data: {
                    type: processNode.spec,
                    value: processNode.codec.encode(data),
                  },
                  inputs: head.process.inputs,
                })
              })
              const res = z.object({ head: z.string(), rebased: z.string() }).parse(await req.json())
              router.push(`${session_id ? `/session/${session_id}` : ''}/report/${res.head}`, undefined, { shallow: true, scroll: false })
            }}
          />
          : <div className="prose max-w-none">Waiting for input</div>}
        {/* <FigureCaption id={head.id} story={story} /> */}
      </div>
      {status ? (
        <div className="alert shadow-lg place-content-start align-middle">
          <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" className="stroke-info shrink-0 w-6 h-6"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"></path></svg>
          <code className="prose max-w-none whitespace-pre-line">{status}</code>
        </div>
      ) : null}
      {outputNode && outputNode.spec === 'Error' && output ? outputNode.view(output) : null}
    </div>
  )
}
