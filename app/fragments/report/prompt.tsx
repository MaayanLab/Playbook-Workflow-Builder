import React from 'react'
import { z } from 'zod'
import { useRouter } from 'next/router'
import { PromptMetaNode } from '@/spec/metanode'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { func_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { Metapath, useMetapathInputs } from '@/app/fragments/metapath'
import type KRG from '@/core/KRG'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function Prompt({ session_id, krg, processNode, output, id, head }: { session_id?: string, krg: KRG, processNode: PromptMetaNode, output: any, id: string, head: Metapath }) {
  const router = useRouter()
  const { data: inputs, error } = useMetapathInputs({ session_id, krg, head })
  const Component = processNode.prompt
  return (
    <div className="collapse collapse-arrow text-black dark:text-white">
      <input type="checkbox" defaultChecked={true} />
      <div className="collapse-title flex flex-row gap-2">
        <Icon icon={processNode.meta.icon || func_icon} className="fill-black dark:fill-white" />
        <h2 className="bp4-heading">{processNode.meta.label || processNode.spec}</h2>
      </div>
      {error ? <div className="alert alert-error prose">{error.toString()}</div> : null}
      <div className="collapse-content">
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
              router.push(`${session_id ? `/session/${session_id}` : ''}/report/${res.head}`, undefined, { shallow: true, scroll: false })
            }}
          />
          : <div>Waiting for input(s)</div>}
      </div>
    </div>
  )
}
