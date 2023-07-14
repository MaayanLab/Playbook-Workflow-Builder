import React from 'react'
import dynamic from 'next/dynamic'
import Head from 'next/head'
import type { GetServerSidePropsContext } from 'next'
import fpprg from '@/app/fpprg'
import krg from '@/app/krg'
import db from '@/app/db'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import { useRouter } from 'next/router'
import { SWRConfig } from 'swr'
import { MetaNode } from '@/spec/metanode'
import fetcher from '@/utils/next-rest-fetcher'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const Graph = dynamic(() => import('@/app/fragments/graph/graph'))
const UserIdentity = dynamic(() => import('@/app/fragments/graph/useridentity'))

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
  const fpl = await fpprg.getFPL(params.graph_id)
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
    const output = await fpprg.getResolved(result.process.id)
    if (output) fallback[`/api/db/process/${result.process.id}/output`] = output.toJSON().data
  }
  const suggestions = await db.objects.suggestion.findMany()
  for (const key in suggestions) {
    const suggestion = suggestions[key]
    let OutputNode = krg.getDataNode(suggestion.output)
    if (OutputNode === undefined) {
      OutputNode = MetaNode(suggestion.output)
        .meta({
          label: `${suggestion.output} (Suggestion)`,
          description: `A data type, suggested as part of ${suggestion.name}`,
          pagerank: -100,
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
      const ProcessNode = MetaNode(suggestion.name)
        .meta({
          label: `${suggestion.name} (Suggestion)`,
          description: suggestion.description,
          pagerank: -100,
        })
        .inputs(suggestion.inputs ?
            dict.init(suggestion.inputs.split(',').map((spec: string, ind: number) =>
              ({ key: ind.toString(), value: krg.getDataNode(spec) })))
            : {} as any)
        .output(OutputNode)
        .prompt((props) => {
          return <div>
            <p>{suggestion.description}</p>
            <p>This was suggested by {suggestion.user ? <UserIdentity user={suggestion.user} /> : <>a playbook partnership user</>}.</p>
          </div>
        })
        .story(props => `It is suggested that "${suggestion.description}" be applied to the inputs: ${suggestion.inputs} to get a ${OutputNode.meta.label}.`)
        .build()
      krg.add(ProcessNode)
    }
  }
  fallback[`/api/suggest`] = suggestions
  return {
    props: {
      fallback,
      extend,
      suggest,
    }
  }
}

export default function App({ fallback, extend, suggest }: { fallback: any, extend: boolean, suggest: boolean }) {
  const router = useRouter()
  const params = ParamType.parse(router.query)
  const graph_id = typeof params !== 'undefined' && 'graph_id' in params && params.graph_id ? params.graph_id : 'start'
  const node_id = typeof params !== 'undefined' && 'node_id' in params && params.node_id ? params.node_id : graph_id
  return (
    <Layout>
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
    </Layout>
  )
}
