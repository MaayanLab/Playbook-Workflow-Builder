import '@blueprintjs/icons/lib/css/blueprint-icons.css'
import '@blueprintjs/select/lib/css/blueprint-select.css'
import '@blueprintjs/core/lib/css/blueprint.css'

import React from 'react'
import dynamic from 'next/dynamic'
import Head from 'next/head'
import type { GetServerSidePropsContext } from 'next'
import fpprg from '@/app/fpprg'
import krg from '@/app/krg'
import db from '@/app/kvdb'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import { useRouter } from 'next/router'
import { SWRConfig } from 'swr'
import type { Metapath } from '@/app/fragments/graph/types'
import { MetaNode } from '@/spec/metanode'

const Header = dynamic(() => import('@/app/fragments/playbook/header'))
const Footer = dynamic(() => import('@/app/fragments/playbook/footer'))
const Graph = dynamic(() => import('@/app/fragments/graph/graph'), { ssr: false })

const ParamType = z.union([
  z.object({ graph_id: z.string(), node_id: z.string() }),
  z.object({ graph_id: z.string() }),
  z.object({}),
  z.undefined(),
])

export async function getServerSideProps(ctx: GetServerSidePropsContext) {
  // we expect a uri of the form /[graph_id][/node_id]
  const params = ParamType.parse(ctx.params)
  const extend = ctx.resolvedUrl.endsWith('extend')
  const suggest = ctx.resolvedUrl.endsWith('suggest')
  if (params === undefined || !('graph_id' in params) || params.graph_id === 'start') {
    return {
      props: { fallback: {}, extend, suggest }
    }
  }
  const fpl = fpprg.getFPL(params.graph_id)
  if (fpl === undefined) {
    return { notFound: true }
  }
  const results = fpl.resolve().map(fpl => fpl.toJSON())
  if ('node_id' in params && !(params.node_id === 'start' || results.some(({ id }) => id === params.node_id))) {
    return { notFound: true }
  }
  const fallback: Record<string, unknown> = {
    [`/api/db/fpl/${params.graph_id}`]: results
  }
  for (const result of results) {
    const output = fpprg.getResolved(result.process.id)
    if (output) fallback[`/api/db/process/${result.process.id}/output`] = output.toJSON().data
  }
  const kvdb: Record<string, string> = {}
  for await (const [ key, value ] of db.iterator() as any) {
    kvdb[key] = value.toString()
    const suggestion = JSON.parse(kvdb[key])
    let OutputNode = krg.getDataNode(suggestion.output)
    if (OutputNode === undefined) {
      OutputNode = MetaNode.createData(suggestion.output)
        .meta({
          label: suggestion.output,
          description: `A data type, suggested as part of ${suggestion.name}`,
        })
        .codec<any>()
        .view((props) => {
          return <div>This data type was suggested as part of {suggestion.name}</div>
        })
        .build()
      krg.add(OutputNode)
    }
    let ProcessNode = krg.getProcessNode(suggestion.name)
    if (ProcessNode === undefined) {
      const ProcessNode = MetaNode.createProcess(suggestion.name)
        .meta({
          label: suggestion.name,
          description: suggestion.description,
        })
        .inputs(suggestion.inputs ?
            dict.init(suggestion.inputs.split(',').map((spec: string, ind: number) =>
              ({ key: ind.toString(), value: krg.getDataNode(spec) })))
            : {} as any)
        .output(OutputNode)
        .prompt((props) => {
          return <div>This was suggested by {suggestion.author_name} &lt;{suggestion.author_email}&gt; ({suggestion.author_org})</div>
        })
        .build()
      krg.add(ProcessNode)
    }
  }
  fallback[`/api/suggest`] = kvdb
  return {
    props: {
      fallback,
      extend,
      suggest,
    }
  }
}

async function fetcher(path: string): Promise<Array<Metapath>> {
  const req = await fetch(path)
  return await req.json()
}

export default function App({ fallback, extend, suggest }: { fallback: any, extend: boolean, suggest: boolean }) {
  const router = useRouter()
  const params = ParamType.parse(router.query)
  const graph_id = typeof params !== 'undefined' && 'graph_id' in params && params.graph_id ? params.graph_id : 'start'
  const node_id = typeof params !== 'undefined' && 'node_id' in params && params.node_id ? params.node_id : graph_id
  return (
    <div className="flex flex-col min-w-screen min-h-screen">
      <Head>
        <title>Playbook</title>
      </Head>

      <Header homepage="/graph" />

      <SWRConfig value={{ fallback, fetcher }}>
        <main className="flex-grow container mx-auto py-4 flex flex-col">
          <Graph
            graph_id={graph_id}
            node_id={node_id}
            extend={extend}
            suggest={suggest}
          />
        </main>
      </SWRConfig>

      <Footer />
    </div>
  )
}
