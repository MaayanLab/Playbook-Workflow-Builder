import React from 'react'
import type { GetServerSidePropsContext } from 'next'
import { useRouter } from 'next/router'
import fpprg from '@/app/fpprg'
import { FPL } from '@/core/FPPRG'
import krg from '@/app/krg'
import { z } from 'zod'
import useSWR, { SWRConfig } from 'swr'

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
    [`/api/db/${id}`]: results
  }
  for (const result of results) {
    const output = fpprg.getResolved(result.process.id)
    if (output)  fallback[`/api/db/${result.id}/output`] = result.process
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
      <Cells id={id} />
    </SWRConfig>
  )
}

function Cells({ id }: { id?: string }) {
  const router = useRouter()
  const { data: metapath, error } = useSWR<Array<Metapath>>(`/api/db/${id}`)
  const head = metapath ? metapath[metapath.length - 1] : undefined
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
          <button
            onClick={async () => {
              const inputs: Record<string, unknown> = {}
              if (head) {
                for (const i in proc.inputs) {
                  inputs[i] = { id: head.process.id }
                }
              }
              const req = await fetch(`/api/db/${id}/extend`, {
                method: 'POST',
                body: JSON.stringify({
                  type: proc.spec,
                  inputs,
                })
              })
              const res = z.string().parse(await req.json())
              router.push(`/persistent/${res}`, undefined, { shallow: true })
            }}
          >{proc.spec}</button>
          <span>&nbsp; =&gt; {proc.output.spec}</span>
        </div>
      ) : null}
    </div>
  )
}

function Cell({ id, head }: { id?: string, head: Metapath }) {
  const router = useRouter()
  const { data: rawOutput, error: outputError } = useSWR(`/api/db/${head.id}/output`)
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
          const req = await fetch(`/api/db/${id}/rebase/${head.process.id}`, {
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
          const res = z.string().parse(await req.json())
          router.push(`/persistent/${res}`, undefined, { shallow: true })
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
