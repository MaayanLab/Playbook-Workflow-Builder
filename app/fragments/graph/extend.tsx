import dynamic from 'next/dynamic'
import { NextRouter, useRouter } from 'next/router'
import type KRG from '@/core/KRG'
import { z } from 'zod'
import { start_icon, rightarrow_icon, func_icon, variable_icon } from '@/icons'
import { MetaNodePromptType, MetaNodeResolveType } from '@/spec/metanode'
import type { Metapath } from '@/app/fragments/graph/types'
import { SuggestionEdges } from '@/app/fragments/graph/suggest'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'

import type CatalogType from '@/app/fragments/graph/catalog'
const Catalog = dynamic(() => import('@/app/fragments/graph/catalog')) as typeof CatalogType
const Icon = dynamic(() => import('@/app/components/icon'))
const Card = dynamic(() => import('@blueprintjs/core').then(({ Card }) => Card))

export default function Extend({ krg, id, head }: { krg: KRG, id: string, head: Metapath }) {
  const router = useRouter()
  const processNode = head ? krg.getProcessNode(head.process.type) : undefined
  return (
    <Catalog<MetaNodePromptType | MetaNodeResolveType & ({} | { onClick: (_: { router: NextRouter, id: string, head: Metapath }) => void })>
      items={[
        ...krg.getNextProcess(processNode ? processNode.output.spec : ''),
        ...SuggestionEdges(processNode ? processNode.output : undefined),
      ]}
      serialize={item => [
        item.spec,
        item.meta.label,
        item.meta.description,
        dict.values(item.meta.tags || {})
          .flatMap(tagGroup => dict.items(tagGroup).filter(({ value }) => value).map(({ key }) => key)),
      ].join(' ')}
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
          {Object.keys(item.inputs).length === 0 ? <Icon icon={start_icon} /> : null}
          {dict.items(item.inputs).map(({ key, value }, i) => (
            <span key={key.toString()}>
              {i > 0 ? <Icon icon={rightarrow_icon} /> : null}
              <Icon icon={array.ensureOne(value).meta.icon || variable_icon} />
            </span>
          ))}
          <Icon icon={rightarrow_icon} />
          <Icon icon={item.meta.icon || func_icon} />
          <Icon icon={rightarrow_icon} />
          <Icon icon={item.output.meta.icon || variable_icon} />
        </div>
        <h5 className="bp4-heading">{item.meta.label || ''}</h5>
        <p className="bp4-text-small">{item.meta.description || ''}</p>
      </Card>
    }</Catalog>
  )
}
