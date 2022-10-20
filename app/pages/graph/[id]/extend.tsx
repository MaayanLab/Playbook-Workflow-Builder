import '@blueprintjs/icons/lib/css/blueprint-icons.css'
import '@blueprintjs/select/lib/css/blueprint-select.css'
import '@blueprintjs/core/lib/css/blueprint.css'

import dynamic from 'next/dynamic'
import Head from 'next/head'
import { useRouter } from 'next/router'
import krg from '@/app/krg'
import { z } from 'zod'

import { start_icon, rightarrow_icon } from '@/icons'
import { MetaNodePromptType, MetaNodeResolveType, MetaNodeType } from '@/spec/metanode'

const Header = dynamic(() => import('@/app/fragments/playbook/header'))
const Footer = dynamic(() => import('@/app/fragments/playbook/footer'))
const Breadcrumbs = dynamic(() => import('@/app/fragments/breadcrumbs'))
import type CatalogType from '@/app/fragments/catalog'
const Catalog = dynamic(() => import('@/app/fragments/catalog')) as typeof CatalogType
const Icon = dynamic(() => import('@/app/components/icon'))
const Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export default function App() {
  const router = useRouter()
  return (
    <div className="flex flex-col min-w-screen min-h-screen">
      <Head>
        <title>Playbook</title>
      </Head>

      <Header />

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
            {
              id: 'extend',
              kind: 'process',
              label: 'Extend',
              color: '#B3CFFF',
              content: '+',
              parents: ['start'],
            }
          ]}
          onclick={(_evt, id) => {
            if (id === 'start') {
              router.push(`/graph`, undefined, { shallow: true })
            }
          }}
        />
      </div>

      <main className="flex-grow m-4">
        <Catalog<MetaNodePromptType | MetaNodeResolveType>
          items={krg.getNextProcess('')}
          // items={Object.values(metagraph.getNeighbors({ graph: ctx.graph, node: ctx.node })) as SCG.MetaEdge[]}
          serialize={(item) => item.spec}
          // serialize={(item: SCG.MetaEdge) => `${ensureCallable((item.meta || {}).name)(ctx)} ${ensureCallable((item.meta || {}).desc)(ctx)} ${ensureCallable(metagraph.nodes[item.input.spec].meta.name)(ctx)} ${ensureCallable(metagraph.nodes[item.output.spec].meta.name)(ctx)}`}
        >{(item) =>
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
                const req = await fetch(`/api/db/fpl`, {
                  method: 'POST',
                  body: JSON.stringify([{
                    type: item.spec,
                    inputs: {}, // TODO
                  }])
                })
                const res = z.string().parse(await req.json())
                router.push(`/graph/${res}`)
              }}
            >
              Select
            </Button>
          </div>
        }</Catalog>
      </main>

      <Footer />
    </div>
  )
}