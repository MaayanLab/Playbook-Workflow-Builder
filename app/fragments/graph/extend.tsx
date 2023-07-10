import React from 'react'
import dynamic from 'next/dynamic'
import { NextRouter, useRouter } from 'next/router'
import type KRG from '@/core/KRG'
import { z } from 'zod'
import { start_icon, rightarrow_icon, func_icon, variable_icon } from '@/icons'
import { ProcessMetaNode } from '@/spec/metanode'
import type { Metapath } from '@/app/fragments/metapath'
import { SuggestionEdges } from '@/app/fragments/graph/suggest'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import Head from 'next/head'

import type CatalogType from '@/app/fragments/graph/catalog'
const Catalog = dynamic(() => import('@/app/fragments/graph/catalog')) as typeof CatalogType
const Icon = dynamic(() => import('@/app/components/icon'))
const Card = dynamic(() => import('@blueprintjs/core').then(({ Card }) => Card))

export default function Extend({ krg, id, head, metapath }: { krg: KRG, id: string, head: Metapath, metapath: Metapath[] }) {
  const router = useRouter()
  const processNode = head ? krg.getProcessNode(head.process.type) : undefined
  const selections = React.useMemo(() => {
    // we'll use leaf nodes of the metapath + the current selected node as the selections
    const selections: Record<string, { process: Metapath["process"], processNode: ProcessMetaNode }> = {}
    ;[...metapath, head].forEach(item => {
      if (item === undefined) return
      // add this to the selections
      selections[item.process.id] = { process: item.process, processNode: krg.getProcessNode(item.process.type) }
      // if a selection previously registered is a parent of this selection, remove it from selections
      dict.values(item.process.inputs).forEach(k => {
        if (k.id in selections) delete selections[k.id]
      })
    })
    return selections
  }, [metapath, head])
  return (
    <>
      <Head>
        <title>Playbook: Extend{processNode?.output.meta.label ? ` from ${processNode?.output.meta.label}` : null}</title>
      </Head>
      <Catalog<ProcessMetaNode & ({}|{ onClick: (_: { router: NextRouter, id: string, head: Metapath }) => void })>
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
                    const relevantSelections = dict.filter(selections, ({ value: selection }) => selection.processNode.output.spec === input.spec)
                    const selection = head.process.id in relevantSelections ? head : array.ensureOne(dict.values(relevantSelections))
                    inputs[arg] = { id: selection.process.id }
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
              {dict.isEmpty(item.inputs) ? <Icon icon={start_icon} /> : null}
              {dict.items(item.inputs).map(({ key, value }, i) => {
                const item_input = array.ensureOne(value)
                return (
                  <span key={key.toString()} className="flex flex-row items-center">
                    {i > 0 ? <Icon title={null} icon={rightarrow_icon} /> : null}
                    {Array.isArray(value) ? '[' : null}
                    <Icon title={item_input.meta.label} icon={item_input.meta.icon || variable_icon} />
                    {Array.isArray(value) ? ']' : null}
                  </span>
                )
              })}
              <Icon title={null} icon={rightarrow_icon} />
              <Icon title={item.meta.label} icon={item.meta.icon || func_icon} />
              <Icon title={null} icon={rightarrow_icon} />
              <Icon title={item.output.meta.label} icon={item.output.meta.icon || variable_icon} />
            </div>
            <h5 className="text-lg font-bold">{item.meta.label || ''}</h5>
            <p className="text-sm">{item.meta.description || ''}</p>
          </Card>
        )
      }}</Catalog>
    </>
  )
}
