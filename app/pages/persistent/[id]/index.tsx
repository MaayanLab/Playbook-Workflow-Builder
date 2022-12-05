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

export default function App({ id, metapath }: { id: string, metapath: Array<UnPromise<ReturnType<FPL['toJSONWithOutput']>>> }) {
  const router = useRouter()
  const head = metapath[metapath.length-1]
  const processNode = krg.getProcessNode(head.process.type)
  if ('prompt' in processNode) {
    const inputs: any = head.process.inputs
    const Prompt = processNode.prompt
    return <Prompt
      inputs={inputs}
      output={head.process.output ? processNode.output.codec.decode(head.process.output.value) : undefined}
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
    />
  }
  return (
    <div>
      {head.process.type} ({head.id})
    </div>
  )
}
