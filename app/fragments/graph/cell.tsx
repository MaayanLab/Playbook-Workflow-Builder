import React from 'react'
import useSWRImmutable from 'swr/immutable'
import { z } from 'zod'
import { useRouter } from 'next/router'
import type { Metapath } from '@/app/fragments/graph/types'
import type KRG from '@/core/KRG'
import { TimeoutError } from '@/spec/error'

export default function Cell({ krg, id, head, autoextend }: { krg: KRG, id: string, head: Metapath, autoextend: boolean }) {
  const router = useRouter()
  const { data: rawOutput, error: outputError } = useSWRImmutable(`/api/db/process/${head.process.id}/output`)
  const processNode = React.useMemo(() => krg.getProcessNode(head.process.type), [head])
  const outputNode = React.useMemo(() => {
    if (!processNode) return
    if (!rawOutput || outputError) return processNode.output
    return krg.getDataNode(rawOutput.type)
  }, [processNode, rawOutput, outputError])
  const inputs = head.process.inputs
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
  const Prompt = 'prompt' in processNode ? processNode.prompt : undefined
  return (
    <div className="flex-grow flex flex-col">
      <div className="mb-4">
        <h2 className="bp4-heading">{processNode.meta.label || processNode.spec}</h2>
        {Prompt ? <Prompt
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
                inputs,
              })
            })
            const res = z.object({ head: z.string(), rebased: z.string() }).parse(await req.json())
            router.push(`/graph/${res.head}${res.head !== res.rebased ? `/node/${res.rebased}` : ''}${autoextend ? '/extend' : ''}`, undefined, { shallow: true })
          }}
        />
        : processNode.meta.description ? <p className="bp4-ui-text">{processNode.meta.description}</p>
        : null}
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
    </div>
  )
}
