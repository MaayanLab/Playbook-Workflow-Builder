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

export default function Extend({ krg, id, head, metapath }: { krg: KRG, id: string, head: Metapath, metapath: Metapath[] }) {
  const router = useRouter()
  const processNode = head ? krg.getProcessNode(head.process.type) : undefined
  const selections = dict.init(metapath.map(item => ({ key: item.process.id, value: { process: item.process, processNode: krg.getProcessNode(item.process.type) } })))
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
    >{item => {
      // determine if multi-inputs are satisfiable
      const disabled = !array.all(
        dict.values(item.inputs).map((value) => {
          if (Array.isArray(value)) {
            return dict.values(selections).filter(selection => selection.processNode.output.spec === value[0].spec).length > 1
          } else {
            return dict.values(selections).filter(selection => selection.processNode.output.spec === value.spec).length >= 1
          }
        })
      )
      return (
        <Card
          key={item.spec}
          interactive={!disabled}
          style={{
            backgroundColor: item.output.meta.color || 'lightgrey',
            opacity: disabled ? 0.75 : 1,
          }}
          onClick={async () => {
            if (disabled) return
            if ('onClick' in item) {
              item.onClick({ router, id, head })
            } else {
              const inputs: Record<string, { id: string }> = {}
              dict.items(item.inputs).forEach(({ key: arg, value: input }) => {
                if (Array.isArray(input)) {
                  dict.values(selections)
                    .filter(selection => selection.processNode.output.spec === input[0].spec)
                    .forEach((selection, i) => {
                      inputs[`${arg}:${i}`] = { id: selection.process.id }
                    })
                } else {
                  dict.values(selections)
                    .filter(selection => selection.processNode.output.spec === input.spec)
                    .forEach(selection => {
                      inputs[arg] = { id: selection.process.id }
                    })
                }
              })
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
              <span key={key.toString()} className="flex flex-row items-center">
                {i > 0 ? <Icon icon={rightarrow_icon} /> : null}
                {Array.isArray(value) ? '[' : null}
                <Icon icon={array.ensureOne(value).meta.icon || variable_icon} />
                {Array.isArray(value) ? ']' : null}
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
      )
    }}</Catalog>
  )
}
