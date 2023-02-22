import React from 'react'
import { z } from 'zod'
import { useRouter } from 'next/router'
import type { Metapath } from '@/app/fragments/graph/types'
import { PromptMetaNode } from '@/spec/metanode'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { func_icon } from '@/icons'
import { useSWRImmutableSticky } from '@/utils/use-sticky'
import dynamic from 'next/dynamic'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function Prompt({ processNode, output, id, head }: { processNode: PromptMetaNode, output: any, id: string, head: Metapath }) {
  const router = useRouter()
  const { data: rawInputs, error: inputsError } = useSWRImmutableSticky(() => !dict.isEmpty(processNode.inputs) ? `/api/db/process/${head.process.id}/inputs` : undefined)
  const { inputs, decodeError } = React.useMemo(() => {
    if (inputsError) return {}
    if (!rawInputs) return { inputs: dict.isEmpty(processNode.inputs) ? {} : undefined }
    try {
      return {
        inputs: dict.init(
          dict.items(processNode.inputs).map(({ key, value }) => {
            if (Array.isArray(value)) {
              return { key, value: rawInputs[key].map(value[0].codec.decode) }
            } else {
              return { key, value: value.codec.decode(rawInputs[key]) }
            }
          })
        ),
      }
    } catch (e: any) {
      console.error(e)
      return { decodeError: process.env.NODE_ENV === 'production' ? 'An unexpected error occurred.' : e.toString() }
    }
  }, [rawInputs, inputsError, processNode])
  const Component = processNode.prompt
  return (
    <div className="p-3">
      <div className="flex flex-row gap-2">
        <Icon icon={processNode.meta.icon || func_icon} />
        <h2 className="bp4-heading">{processNode.meta.label || processNode.spec}</h2>
      </div>
      {decodeError ? <div className="alert alert-error">{decodeError.toString()}</div> : null}
      {inputs !== undefined && array.intersection(dict.keys(processNode.inputs), dict.keys(inputs)).length === dict.keys(processNode.inputs).length ?
        <Component
          inputs={inputs}
          output={output}
          submit={async (output) => {
            const req = await fetch(`/api/db/fpl/${id}/rebase/${head.process.id}`, {
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
            router.push(`/report/${res.head}`, undefined, { shallow: true, scroll: false })
          }}
        />
        : <div>Waiting for input(s)</div>}
    </div>
  )
}
