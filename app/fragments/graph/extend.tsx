import dynamic from 'next/dynamic'
import { NextRouter, useRouter } from 'next/router'
import krg from '@/app/krg'
import { z } from 'zod'
import { start_icon, rightarrow_icon } from '@/icons'
import { MetaNodePromptType, MetaNodeResolveType } from '@/spec/metanode'
import type { Metapath } from '@/app/fragments/graph/types'
import { SuggestionEdges } from '@/app/fragments/graph/suggest'

import type CatalogType from '@/app/fragments/graph/catalog'
const Catalog = dynamic(() => import('@/app/fragments/graph/catalog')) as typeof CatalogType
const Icon = dynamic(() => import('@/app/components/icon'))
const Card = dynamic(() => import('@blueprintjs/core').then(({ Card }) => Card))

export default function Extend({ id, head }: { id: string, head: Metapath }) {
  const router = useRouter()
  const processNode = head ? krg.getProcessNode(head.process.type) : undefined
  return (
    <Catalog<MetaNodePromptType | MetaNodeResolveType & ({} | { onClick: (_: { router: NextRouter, id: string, head: Metapath }) => void })>
      items={[
        ...krg.getNextProcess(processNode ? processNode.output.spec : ''),
        ...SuggestionEdges(processNode ? processNode.output : undefined),
      ]}
      // items={Object.values(metagraph.getNeighbors({ graph: ctx.graph, node: ctx.node })) as SCG.MetaEdge[]}
      serialize={item => item.spec}
      // serialize={(item: SCG.MetaEdge) => `${ensureCallable((item.meta || {}).name)(ctx)} ${ensureCallable((item.meta || {}).desc)(ctx)} ${ensureCallable(metagraph.nodes[item.input.spec].meta.name)(ctx)} ${ensureCallable(metagraph.nodes[item.output.spec].meta.name)(ctx)}`}
    >{item =>
      <Card
        key={item.spec}
        interactive={true}
        style={{
          backgroundColor: item.output.meta.color || 'lightgrey',
        }}
        onClick={async () => {
          if ('onClick' in item) {
            item.onClick({ router, id, head })
          } else {
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
          }
        }}
      >
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
        <h5 className="bp4-heading">{item.meta.label || ''}</h5>
        <p className="bp4-text-small">{item.meta.description || ''}</p>
      </Card>
    }</Catalog>
  )
}
