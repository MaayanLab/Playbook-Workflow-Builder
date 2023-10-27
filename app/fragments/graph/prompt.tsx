import React from 'react'
import { z } from 'zod'
import { useRouter } from 'next/router'
import { DataMetaNode, PromptMetaNode } from '@/spec/metanode'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import type KRG from '@/core/KRG'
import { Metapath, useMetapathInputs } from '@/app/fragments/metapath'
import { useStory } from '@/app/fragments/story'
import { TimeoutError } from '@/spec/error'

export default function Prompt({ session_id, krg, processNode, outputNode, output, id, head, autoextend }: { session_id?: string, krg: KRG, processNode: PromptMetaNode, outputNode: DataMetaNode, output: any, id: string, head: Metapath, autoextend: boolean }) {
  const router = useRouter()
  const { data: inputs, error } = useMetapathInputs({ session_id, krg, head })
  const story = useStory()
  const [storyText, storyCitations] = React.useMemo(() => story.split('\n\n'), [story])
  const Component = processNode.prompt
  return (
    <div className="flex-grow flex flex-col">
      <div className="mb-4">
        <h2 className="bp5-heading">{processNode.meta.label || processNode.spec}</h2>
        <p className="prose">{storyText}</p>
        <p className="prose text-sm">{storyCitations}</p>
      </div>
      {error && !(error instanceof TimeoutError) ? <div className="alert alert-error prose">{error.toString()}</div> : null}
      {outputNode && outputNode.spec === 'Error' && output ? outputNode.view(output) : null}
      {inputs !== undefined && array.intersection(dict.keys(processNode.inputs), dict.keys(inputs)).length === dict.keys(processNode.inputs).length ?
        <Component
          session_id={session_id}
          data={head.process.data?.value}
          inputs={inputs}
          output={output}
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
            router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${res.head}${res.head !== res.rebased ? `/node/${res.rebased}` : ''}${autoextend ? '/extend' : ''}`, undefined, { shallow: true })
          }}
        />
        : <div>Waiting for input</div>}
    </div>
  )
}
