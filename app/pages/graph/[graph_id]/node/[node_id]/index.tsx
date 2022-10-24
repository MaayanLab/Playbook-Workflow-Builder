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
import { rightarrow_icon, start_icon } from '@/icons'
import { MetaNodePromptType, MetaNodeResolveType } from '@/spec/metanode'

const Header = dynamic(() => import('@/app/fragments/playbook/header'))
const Footer = dynamic(() => import('@/app/fragments/playbook/footer'))
const Home = dynamic(() => import('@/app/fragments/playbook/home'))
const Breadcrumbs = dynamic(() => import('@/app/fragments/breadcrumbs'))

import type CatalogType from '@/app/fragments/catalog'
const Catalog = dynamic(() => import('@/app/fragments/catalog')) as typeof CatalogType
const Icon = dynamic(() => import('@/app/components/icon'))
const Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

type Metapath = ReturnType<FPL['toJSON']>

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
  if (params === undefined || !('graph_id' in params) || params.graph_id === 'start') {
    return {
      props: { fallback: {}, extend }
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
  return {
    props: {
      fallback,
      extend,
    }
  }
}

async function fetcher(path: string): Promise<Array<Metapath>> {
  const req = await fetch(path)
  return await req.json()
}

export default function App({ fallback, extend }: { fallback: any, extend: boolean }) {
  const router = useRouter()
  const params = ParamType.parse(router.query)
  const graph_id = typeof params !== 'undefined' && 'graph_id' in params && params.graph_id ? params.graph_id : 'start'
  const node_id = typeof params !== 'undefined' && 'node_id' in params && params.node_id ? params.node_id : graph_id
  return (
    <div className="flex flex-col min-w-screen min-h-screen">
      <Head>
        <title>Playbook</title>
      </Head>

      <Header />

      <SWRConfig value={{ fallback, fetcher }}>
        <main className="flex-grow m-4 flex flex-col">
          <Graph
            graph_id={graph_id}
            node_id={node_id}
            extend={extend}
          />
        </main>
      </SWRConfig>

      <Footer />
    </div>
  )
}

function Graph({ graph_id, node_id, extend }: { graph_id: string, node_id: string, extend: boolean }) {
  const router = useRouter()
  const { data: metapath_, error } = useSWRImmutable<Array<Metapath>>(() => graph_id !== 'start' ? `/api/db/fpl/${graph_id}` : undefined)
  const metapath = metapath_ || []
  const head = metapath.filter(({ id }) => id === node_id)[0]
  return (
    <>
      <div className="flex w-auto h-40">
        <Breadcrumbs
          graph={[
            {
              id: 'start',
              kind: 'data',
              label: 'Start',
              color: node_id === 'start' ? '#B3CFFF' : 'lightgrey',
              icon: [start_icon],
              parents: [],
            },
            ...(metapath||[]).flatMap((head, i) => [
              {
                id: head.process.id,
                kind: 'process' as 'process',
                label: head.process.type,
                color: head.id === node_id ? '#B3CFFF' : 'lightgrey',
                icon: [],
                parents: Object.keys(head.process.inputs).length === 0 ? ['start'] : Object.values(head.process.inputs).map(({ id }) => `${id}:output`),
              },
              {
                id: `${head.process.id}:output`,
                kind: 'data' as 'data',
                label: '',
                color: head.id === node_id ? '#B3CFFF' : 'lightgrey',
                icon: [],
                parents: [head.process.id],
              },
            ]),
            {
              id: 'extend',
              kind: 'process' as 'process',
              label: 'Extend',
              color: extend ? '#B3CFFF' : 'lightgrey',
              content: '+',
              parents: [head ? `${head.process.id}:output` : 'start'],
            },
          ]}
          onclick={(_evt, id) => {
            if (id === 'extend') {
              router.push(`/graph/${graph_id}${graph_id !== node_id ? `/node/${node_id}` : ''}/extend`, undefined, { shallow: true })
            } else {
              const focus_node_process_id = id.split(':')[0]
              const focus_node_id = id === 'start' ? 'start' : metapath.filter((node) => node.process.id === focus_node_process_id)[0].id
              router.push(`/graph/${graph_id}${graph_id !== focus_node_id ? `/node/${focus_node_id}` : ''}`, undefined, { shallow: true })
            }
          }}
        />
      </div>
      <main className="flex-grow flex flex-col">
        {error ? <div>{error}</div> : null}
        {extend ?
          <Extend id={graph_id} head={head} />
          : node_id === 'start' ?
            <Home />
            : head ?
              <Cell id={graph_id} head={head} />
              : null}
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
    <div className="flex-grow flex flex-col">
      <div className="flex-grow flex-shrink">
        <h2 className="text-lg mb-1">Process ({processNode.spec})</h2>
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
            router.push(`/graph/${res}`, undefined, { shallow: true })
          }}
        /> : null}
      </div>
      <div className="flex-grow">
        {outputNode ? (
          <>
            <h2 className="text-lg mb-1">View ({outputNode.spec})</h2>
            {View && output ? View(output) : null}
          </>
        ) : (
          <div>Loading...</div>
        )}
      </div>
    </div>
  )
}

function Extend({ id, head }: { id: string, head: Metapath }) {
  const router = useRouter()
  const processNode = head ? krg.getProcessNode(head.process.type) : undefined
  return (
    <Catalog<MetaNodePromptType | MetaNodeResolveType>
      items={krg.getNextProcess(processNode ? processNode.output.spec : '')}
      // items={Object.values(metagraph.getNeighbors({ graph: ctx.graph, node: ctx.node })) as SCG.MetaEdge[]}
      serialize={item => item.spec}
      // serialize={(item: SCG.MetaEdge) => `${ensureCallable((item.meta || {}).name)(ctx)} ${ensureCallable((item.meta || {}).desc)(ctx)} ${ensureCallable(metagraph.nodes[item.input.spec].meta.name)(ctx)} ${ensureCallable(metagraph.nodes[item.output.spec].meta.name)(ctx)}`}
    >{item =>
      <div key={item.spec} className="rounded-lg p-1 m-1" style={{ backgroundColor: item.output.meta.color || 'lightgrey' }}>
        <div className="flex flex-row">
          {Object.keys(item.inputs).length === 0 ? (
            <>
              <Icon icon={start_icon} />
              <Icon icon={rightarrow_icon} />
            </>
          ) : Object.values(item.inputs).map(input => (
            <>
              <Icon icon={input.meta.icon} />
              {'icon' in item.meta ? (
                <>
                  <Icon icon={rightarrow_icon} />
                  <Icon icon={item.meta.icon} />
                </>
              ) : null}
            </>
          ))}
        </div>
        <h5 className="card-title mt-4">{item.meta.label || ''}</h5>
        <p className="card-text">{item.meta.description || ''}</p>
        <Button
          onClick={async () => {
            const inputs: Record<string, { id: string }> = {}
            if (head) {
              for (const arg in item.inputs) {
                inputs[arg] = { id: head.process.id }
              }
            }
            const req = await fetch(`/api/db/fpl/${id}/extend`, {
              method: 'POST',
              body: JSON.stringify({
                type: item.spec,
                inputs,
              })
            })
            const res = z.string().parse(await req.json())
            router.push(`/graph/${res}`)
          }}
        >
          Select
        </Button>
      </div>
    }</Catalog>
  )
}