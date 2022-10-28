import useSWRImmutable from 'swr/immutable'
import { z } from 'zod'
import krg from '@/app/krg'
import { useRouter } from 'next/router'
import type { Metapath } from '@/app/fragments/graph/types'

export default function Cell({ id, head }: { id: string, head: Metapath }) {
  const router = useRouter()
  const { data: rawOutput, error: outputError } = useSWRImmutable(`/api/db/process/${head.process.id}/output`)
  const processNode = krg.getProcessNode(head.process.type)
  const inputs: any = head.process.inputs
  const outputNode = rawOutput && !outputError ? krg.getDataNode(rawOutput.type) : undefined
  const output = outputNode ? outputNode.codec.decode(rawOutput.value) : undefined
  const View = outputNode ? outputNode.view : undefined
  const Prompt = 'prompt' in processNode ? processNode.prompt : undefined
  return (
    <div className="flex-grow flex flex-col">
      <div className="flex-grow flex-shrink">
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
            router.push(`/graph/${res.head}${res.head !== res.rebased ? `/node/${res.rebased}` : ''}`, undefined, { shallow: true })
          }}
        />
        : processNode.meta.description ? <p className="bp4-ui-text">{processNode.meta.description}</p>
        : null}
      </div>
      <div className="flex-grow">
        {outputNode ? (
          <>
            <h2 className="bp4-heading">{outputNode.meta.label || outputNode.spec}</h2>
            {View && output ? View(output) : null}
          </>
        ) : (
          <div>Loading...</div>
        )}
      </div>
    </div>
  )
}
