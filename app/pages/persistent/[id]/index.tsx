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
      metapath: fpl.resolve().map(fpl => fpl.toJSON()),
    }
  }
}

export default function App({ id, metapath }: { id: string, metapath: Array<ReturnType<FPL['toJSON']>> }) {
  const router = useRouter()
  const current = metapath[metapath.length-1]
  const processNode = krg.getProcessNode(current.process['@type'])
  if ('prompt' in processNode) {
    const inputs: any = current.process.inputs
    const Prompt = processNode.prompt
    return <Prompt
      inputs={inputs}
      submit={async (output) => {
        const req = await fetch(`/api/db/${id}/resolve`, {
          method: 'POST',
          body: JSON.stringify(processNode.output.codec.encode(output))
        })
        await req.json()
      }}
    />
  }
  return (
    <div>
      {current.process['@type']} ({current['@id']})
    </div>
  )
}
