import React from 'react'
import dynamic from 'next/dynamic'
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
import { Metapath, MetapathProvider, useFPL, useMetapathOutput } from '@/app/fragments/metapath'
import { StoryProvider } from '@/app/fragments/story'
import SafeRender from '@/utils/saferender'
import KRG from '@/core/KRG'
import useKRG from '@/app/fragments/graph/krg'

const UserIdentity = dynamic(() => import('@/app/fragments/graph/useridentity'))

const ParamType = z.union([
  z.object({ session_id: z.string().optional(), graph_id: z.string(), node_id: z.string() }),
  z.object({ session_id: z.string().optional(), graph_id: z.string() }),
  z.object({ session_id: z.string().optional() }),
  z.undefined(),
])

export async function getServerSideProps(ctx: GetServerSidePropsContext) {
  // we expect a uri of the form /[graph_id][/node_id]
  const params = ParamType.parse(ctx.params)
  if (params === undefined || !('graph_id' in params) || params.graph_id === 'start' || params.session_id !== undefined) {
    return {
      props: { fallback: {} }
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
              ({ key: ind.toString(), value: krg.getDataNode(spec) }))
              .filter(({ key, value }) => !!value))
            : {} as any)
        .output(OutputNode)
        .prompt((props) => {
          return <div>
            <p>{suggestion.description}</p>
            <p>This was suggested by {suggestion.user ? <UserIdentity user={suggestion.user} /> : <>a playbook workflow builder user</>}.</p>
          </div>
        })
        .story(props => ({ abstract: `It is suggested that "${suggestion.description}" be applied to the inputs: ${suggestion.inputs} to get a ${OutputNode.meta.label}.` }))
        .build()
      krg.add(ProcessNode)
    }
  }
  fallback[`/api/suggest`] = suggestions
  return {
    props: {
      fallback: JSON.parse(JSON.stringify(fallback)),
    }
  }
}

function Cell({ krg, head }: { krg: KRG, head: Metapath }) {
  const { data: { output, outputNode }, status, error: outputError, mutate } = useMetapathOutput({ krg, head })
  return (

    !outputNode ? <div>Loading...</div>
      : !outputNode?.view || output === undefined ? <div className="prose">Loading...</div>
      : output === null ? <div className="prose max-w-none">Waiting for input</div>
      : <SafeRender component={outputNode.view} props={output} />
  )
}

function Main({ session_id, graph_id, node_id }: { session_id?: string, graph_id: string, node_id: string }) {
  const krg = useKRG({ session_id })
  const { data: metapath = [] } = useFPL(graph_id)
  const head = React.useMemo(() => metapath.filter(({ id }) => id === node_id)[0], [metapath, node_id])
  return (
    <main className="flex-grow flex flex-col p-2">
      <StoryProvider krg={krg} metapath={metapath}>
        {head ? <Cell krg={krg} head={head} /> : null}
      </StoryProvider>
    </main>
  )
}

export default function App({ fallback }: { fallback: any }) {
  const router = useRouter()
  const params = ParamType.parse(router.query)
  const graph_id = typeof params !== 'undefined' && 'graph_id' in params && params.graph_id ? params.graph_id : 'start'
  const node_id = typeof params !== 'undefined' && 'node_id' in params && params.node_id ? params.node_id : graph_id
  return (
    <SWRConfig value={{ fallback, fetcher }}>
      <MetapathProvider session_id={params?.session_id}>
        <Main session_id={params?.session_id} graph_id={graph_id} node_id={node_id} />
      </MetapathProvider>
    </SWRConfig>
  )
}
