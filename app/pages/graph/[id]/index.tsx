import '@blueprintjs/icons/lib/css/blueprint-icons.css'
import '@blueprintjs/select/lib/css/blueprint-select.css'
import '@blueprintjs/core/lib/css/blueprint.css'

import React from 'react'
import dynamic from 'next/dynamic'
import Head from 'next/head'
import type { GetServerSidePropsContext } from 'next'
import { useRouter } from 'next/router'
import fpprg from '@/app/fpprg'
import { FPL } from '@/core/FPPRG'
import krg from '@/app/krg'
import { z } from 'zod'
import useSWRImmutable from 'swr/immutable'
import { SWRConfig } from 'swr'
import { start_icon } from '@/icons'

const Header = dynamic(() => import('@/app/fragments/header'))
const Footer = dynamic(() => import('@/app/fragments/footer'))
const Breadcrumbs = dynamic(() => import('@/app/fragments/breadcrumbs'))

type Metapath = ReturnType<FPL['toJSON']>

const ParamType = z.object({
  id: z.string(),
})

const QueryType = z.object({
  id: z.union([z.string(), z.undefined()]),
})

export async function getServerSideProps(ctx: GetServerSidePropsContext) {
  const { id } = ParamType.parse(ctx.params)
  const fpl = fpprg.getFPL(id)
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
    <div className="flex flex-col min-w-screen min-h-screen">
      <Head>
        <title>Playbook</title>
      </Head>

      <Header />

      <SWRConfig value={{ fallback, fetcher }}>
        <main className="flex-grow m-4">
          <Graph id={id} />
        </main>
      </SWRConfig>

      <Footer />
    </div>
  )
}

function Graph({ id }: { id?: string }) {
  const router = useRouter()
  const { data: metapath, error } = useSWRImmutable<Array<Metapath>>(id ? `/api/db/fpl/${id}` : undefined)
  if (!metapath) return null
  const head = metapath[metapath.length-1]
  return (
    <>
      <div className="flex w-auto h-40">
        <Breadcrumbs
          graph={[
            {
              id: 'start',
              kind: 'data',
              label: 'Start',
              color: 'lightgrey',
              icon: [start_icon],
              parents: [],
            },
            ...metapath.flatMap((head, i) => [
              {
                id: head.process.id,
                kind: 'process' as 'process',
                label: head.process.type,
                color: i === metapath.length-1 ? '#B3CFFF' : 'lightgrey',
                icon: [],
                parents: Object.keys(head.process.inputs).length === 0 ? ['start'] : Object.values(head.process.inputs).map(({ id }) => id),
              },
              {
                id: `${head.process.id}:output`,
                kind: 'data' as 'data',
                label: '',
                color: i === metapath.length-1 ? '#B3CFFF' : 'lightgrey',
                icon: [],
                parents: [head.process.id],
              },
            ]),
            {
              id: 'extend',
              kind: 'process',
              label: 'Extend',
              color: 'lightgrey',
              content: '+',
              parents: [`${head.process.id}:output`],
            }
          ]}
          onclick={(_evt, id) => {
            if (id === 'start') {
              router.push(`/graph`, undefined, { shallow: true })
            }
          }}
        />
      </div>
      <main className="m-4">
        {error ? <div>{error}</div> : null}
        {id && head ? <Cell id={id} head={head} /> : null}
      </main>
    </>
  )
}

function Cell({ id, head }: { id: string, head: Metapath }) {
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
