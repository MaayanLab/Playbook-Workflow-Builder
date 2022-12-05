import React from 'react'
import { GetServerSideProps } from 'next'
import { useRouter } from 'next/router'
import fpprg from '@/app/fpprg'
import { FPL } from '@/core/FPPRG'
import krg from '@/app/krg'
import { z } from 'zod'

const ParamType = z.object({
  id: z.string(),
})

export const getServerSideProps: GetServerSideProps = async (ctx) => {
  ctx.res.setHeader(
    'Cache-Control',
    'public, s-maxage=10, stale-while-revalidate=59'
  )
  const { id } = ParamType.parse(ctx.params)
  const fpl = fpprg.getFPL(id)
  if (fpl === undefined) {
    return { notFound: true }
  }
  return {
    props: {
      id,
      metapath: await Promise.all(fpl.resolve().map(fpl => fpl.toJSONWithOutput())),
    }
  }
}

type UnPromise<T> = T extends Promise<infer t> ? t : never

function Cell({ id, head }: { id: string, head: UnPromise<ReturnType<FPL['toJSONWithOutput']>> }) {
  const router = useRouter()
  const processNode = krg.getProcessNode(head.process.type)
  const inputs: any = head.process.inputs
  const outputNode = head.process.output ? krg.getDataNode(head.process.output.type) : undefined
  const output = head.process.output && outputNode ? outputNode.codec.decode(head.process.output.value) : undefined
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
          router.push(`/persistent/${res}`)
        }}
      /> : null}
      {outputNode ? <h2>View ({outputNode.spec})</h2> : null}
      {View && output ? View(output) : null}
    </div>
  )
}

export default function App({ id, metapath }: { id: string, metapath: Array<UnPromise<ReturnType<FPL['toJSONWithOutput']>>> }) {
  const router = useRouter()
  const head = metapath[metapath.length-1]
  return (
    <div>
      {metapath.map(head => <Cell key={head.id} id={id} head={head} />)}
      <h2>Actions</h2>
      {krg.getNextProcess(head.process.output ? head.process.output.type : '').map(proc =>
        <div key={proc.spec}>
          {Object.keys(proc.inputs).length > 0 ?
            <span>{Object.values(proc.inputs).map((i) => i.spec).join(', ')} =&gt;&nbsp;</span>
            : null}
          <button
            onClick={async () => {
              const inputs: Record<string, unknown> = {}
              for (const i in proc.inputs) {
                inputs[i] = { id: head.process.id }
              }
              const req = await fetch(`/api/db/${id}/extend`, {
                method: 'POST',
                body: JSON.stringify({
                  type: proc.spec,
                  inputs,
                })
              })
              const res = z.string().parse(await req.json())
              router.push(`/persistent/${res}`)
            }}
          >{proc.spec}</button>
          <span>&nbsp; =&gt; {proc.output.spec}</span>
        </div>
      )}
    </div>
  )
}
