import React from 'react'
import type { GetServerSidePropsContext } from 'next'
import { useRouter } from 'next/router'
import fpprg from '@/app/fpprg'
import { FPL } from '@/core/FPPRG'
import krg from '@/app/krg'
import { z } from 'zod'
import useSWRImmutable from 'swr/immutable'
import { SWRConfig } from 'swr'
import dynamic from 'next/dynamic'
import type KRG from '@/core/KRG'

const Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

type Metapath = ReturnType<FPL['toJSON']>

const QueryType = z.object({
  id: z.union([z.string(), z.undefined()]),
})

export async function getServerSideProps(ctx: GetServerSidePropsContext) {
  const { id } = QueryType.parse(ctx.params)
  const fpl = id ? fpprg.getFPL(id) : undefined
  if (fpl === undefined) {
    return { notFound: true }
  }
  const results = fpl.resolve().map(fpl => fpl.toJSON())
  const fallback: Record<string, unknown> = {
    [`/api/db/fpl/${id}`]: results
  }
  for (const result of results) {
    const output = fpprg.getResolved(result.process.id)
    if (output) fallback[`/api/db/process/${result.process.id}/output`] = output.toJSON().data
  }
  return {
    props: {
      fallback,
    }
  }
}

async function fetcher(path: string): Promise<Array<Metapath>> {
  const req = await fetch(path)
  return await req.json()
}


export default function App({ fallback }: { fallback: any }) {
  const router = useRouter()
  const { id } = QueryType.parse(router.query)
  return (
    <SWRConfig value={{ fallback, fetcher }}>
      <Cells krg={krg} id={id} />
    </SWRConfig>
  )
}

function Cells({ krg, id }: { krg: KRG, id?: string }) {
  const router = useRouter()
  const { data: metapath_, error } = useSWRImmutable<Array<Metapath>>(id ? `/api/db/fpl/${id}` : undefined)
  const metapath = metapath_ || []
  const head = metapath[metapath.length - 1]
  const processNode = head ? krg.getProcessNode(head.process.type) : undefined
  return (
    <div>
      {error ? <div>{error}</div> : null}
      {metapath ? metapath.map((head, index) =>
        <Cell key={index} id={id} head={head} />
      ) : null}
      <h2>Actions</h2>
      {processNode ? krg.getNextProcess(processNode.output.spec).map(proc =>
        <div key={proc.spec}>
          {Object.keys(proc.inputs).length > 0 ?
            <span>{Object.values(proc.inputs).map((i) => i.spec).join(', ')} =&gt;&nbsp;</span>
            : null}
          <Button
            onClick={async () => {
              const inputs: Record<string, unknown> = {}
              if (head) {
                for (const i in proc.inputs) {
                  inputs[i] = { id: head.process.id }
                }
              }
              const req = await fetch(`/api/db/fpl/${id}/extend`, {
                method: 'POST',
                body: JSON.stringify({
                  type: proc.spec,
                  inputs,
                })
              })
              const res = z.string().parse(await req.json())
              router.push(`/report/${res}`, undefined, { shallow: true })
            }}
          >{proc.spec}</Button>
          <span>&nbsp; =&gt; {proc.output.spec}</span>
        </div>
      ) : null}
    </div>
  )
}

function Cell({ id, head }: { id?: string, head: Metapath }) {
  const router = useRouter()
  const { data: rawOutput, error: outputError } = useSWRImmutable(`/api/db/process/${head.process.id}/output`)
  const processNode = krg.getProcessNode(head.process.type)
  const inputs: any = head.process.inputs
  const outputNode = rawOutput && !outputError ? krg.getDataNode(rawOutput.type) : undefined
  const output = outputNode ? outputNode.codec.decode(rawOutput.value) : undefined
  const View = outputNode ? outputNode.view : undefined
  const Prompt = 'prompt' in processNode ? processNode.prompt : undefined
  return (
    <div>
      <h2>Process ({processNode.spec})</h2>
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
          router.push(`/report/${res.head}`, undefined, { shallow: true })
        }}
      /> : null}
      {outputNode ? (
        <>
          <h2>View ({outputNode.spec})</h2>
          {View && output ? View(output) : null}
        </>
      ) : (
        <div>Loading...</div>
      )}
    </div>
  )
}
