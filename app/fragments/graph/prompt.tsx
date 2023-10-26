import React from 'react'
import { z } from 'zod'
import { useRouter } from 'next/router'
import { PromptMetaNode } from '@/spec/metanode'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import type KRG from '@/core/KRG'
import { Metapath, useMetapathInputs } from '@/app/fragments/metapath'
import { useStory } from '@/app/fragments/story'

export default function Prompt({ session_id, krg, processNode, output, id, head, autoextend }: { session_id?: string, krg: KRG, processNode: PromptMetaNode, output: any, id: string, head: Metapath, autoextend: boolean }) {
  const router = useRouter()
  const { data: inputs, error } = useMetapathInputs({ session_id, krg, head })
  const story = useStory()
  const [storyText, storyCitations] = React.useMemo(() => story.split('\n\n'), [story])
  const Component = processNode.prompt
  return (
    <div className="flex-grow flex flex-col">
      <div className="mb-4">
        <h2 className="bp4-heading">{processNode.meta.label || processNode.spec}</h2>
        <p className="prose">{storyText}</p>
        <p className="prose text-sm">{storyCitations}</p>
      </div>
      {error ? <div className="alert alert-error prose">{error.toString()}</div> : null}
      {inputs !== undefined && array.intersection(dict.keys(processNode.inputs), dict.keys(inputs)).length === dict.keys(processNode.inputs).length ?
        <Component
          session_id={session_id}
          inputs={inputs}
          output={output}
          submit={async (output) => {
            const req = await fetch(`${session_id ? `/api/socket/${session_id}` : ''}/api/db/fpl/${id}/rebase/${head.process.id}`, {
              headers: {
                'Content-Type': 'application/json',
              },
              method: 'POST',
              body: JSON.stringify({
                type: head.process.type,
                data: {
                  type: processNode.output.spec,
                  value: processNode.output.codec.encode(output),
                },
                inputs: head.process.inputs,
              })
            })
            const res = z.object({ head: z.string(), rebased: z.string() }).parse(await req.json())
            router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${res.head}${res.head !== res.rebased ? `/node/${res.rebased}` : ''}${autoextend ? '/extend' : ''}`, undefined, { shallow: true })
          }}
        />
        : <div>Waiting for input(s)</div>}
    </div>
  )
}
