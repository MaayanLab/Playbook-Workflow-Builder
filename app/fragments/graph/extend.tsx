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
import pathWeights from '@/app/public/weights.json'
import Head from 'next/head'
import classNames from 'classnames'
import pluralize from 'pluralize'

import type CatalogType from '@/app/fragments/graph/catalog'
const Catalog = dynamic(() => import('@/app/fragments/graph/catalog')) as typeof CatalogType
const Icon = dynamic(() => import('@/app/components/icon'))
const Card = dynamic(() => import('@blueprintjs/core').then(({ Card }) => Card))

export default function Extend({ session_id, krg, id, heads, metapath }: { session_id?: string, krg: KRG, id: string, heads: Metapath[], metapath: Metapath[] }) {
  const router = useRouter()
  const processNode = heads[0] ? krg.getProcessNode(heads[0].process.type) : undefined
  const { items, selections } = React.useMemo(() => {
    // we'll use leaf nodes of the metapath + the current selected node as the selections
    const selections: Record<string, { process: Metapath["process"], processNode: ProcessMetaNode }> = {}
    if (heads.length <= 1) {
      ;[...metapath, heads[0]].forEach(item => {
        if (item === undefined) return
        // add this to the selections
        selections[item.process.id] = { process: item.process, processNode: krg.getProcessNode(item.process.type) }
        // if a selection previously registered is a parent of this selection, remove it from selections
        dict.values(item.process.inputs).forEach(k => {
          if (k.id in selections) delete selections[k.id]
        })
      })
      const items = [
        ...krg.getNextProcess(processNode ? processNode.output.spec : '')
          .filter(proc => proc.meta.hidden !== true)
          .map(proc => ({ ...proc, meta: { ...proc.meta, tags: { ...(proc.meta.tags??{}), External: { [proc.meta.external ? 'True': 'False']: 1 } } } })),
        ...SuggestionEdges(processNode ? processNode.output : undefined),
      ]
      return { items, selections }
    } else {
      heads.forEach(item => {
        if (item === undefined) return
        selections[item.process.id] = { process: item.process, processNode: krg.getProcessNode(item.process.type) }
      })
      const items = krg.getNextProcess(processNode ? processNode.output.spec : '')
        .filter(proc => proc.meta.hidden !== true && dict.values(proc.inputs).some(value => Array.isArray(value)))
        .map(proc => ({ ...proc, meta: { ...proc.meta, tags: { ...(proc.meta.tags??{}), External: { [proc.meta.external ? 'True': 'False']: 1 } } } }))
      return {
        selections,
        items,
      }
    }
  }, [metapath, heads])
  const weights = React.useMemo(() => {
    const key = ['Start', ...metapath.map(p => p.process.type)].slice(-2).join(' ')
    const weights = pathWeights[key as keyof typeof pathWeights] || {}
    return weights
  }, [metapath])
  return (
    <>
      <Head>
        <title>Playbook: Extend{processNode?.output.meta.label ? ` from ${processNode?.output.meta.label}` : null}</title>
      </Head>
      <Catalog<ProcessMetaNode & ({}|{ onClick: (_: { router: NextRouter, id: string, head: Metapath }) => void })>
        items={items}
        serialize={item => [
          item.spec,
          item.meta.label,
          item.meta.description,
          dict.values(item.meta.tags || {})
            .flatMap(tagGroup => dict.items(tagGroup).filter(({ value }) => value).map(({ key }) => key)),
        ].join(' ')}
        weights={weights}
      >{item => {
        // determine if multi-inputs are satisfiable, and record a reason if not
        const inputUnsatisfaction = dict.items(item.inputs).map(({ key, value }) => {
          if (Array.isArray(value)) {
            return {
              key,
              value: dict.values(selections).filter(selection => selection.processNode.output.spec === value[0].spec).length <= 1 ? `Multiple ${pluralize(value[0].meta.label)} required` : undefined,
            }
          } else {
            return {
              key,
              value: dict.values(selections).filter(selection => selection.processNode.output.spec === value.spec).length < 1 ? `${value.meta.label} required` : undefined,
            }
          }
        }).filter(({ key, value }) => value !== undefined)
        // if we have no input unsatisfiability reasons, we should be enabled
        const disabled = inputUnsatisfaction.length !== 0
        return (
          <div
            key={item.spec}
            className="tooltip"
            data-tip={disabled ? `${inputUnsatisfaction.map(({ key, value }) => value).join('. ')}` : undefined}
          >
            <Card
              interactive={!disabled}
              style={{
                backgroundColor: item.output.meta.color || 'lightgrey',
                opacity: disabled ? 0.75 : 1,
              }}
              className={classNames('text-left', { 'pointer-events-none': disabled })}
              onClick={async () => {
                if (disabled) return
                if ('onClick' in item) {
                  item.onClick({ router, id, head: heads[0] })
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
                      const selection = heads[0].process.id in relevantSelections ? heads[0] : array.ensureOne(dict.values(relevantSelections))
                      inputs[arg] = { id: selection.process.id }
                    }
                  })
                  const req = await fetch(`${session_id ? `/api/socket/${session_id}` : ''}/api/db/fpl/${id}/extend`, {
                    headers: {
                      'Content-Type': 'application/json',
                    },
                    method: 'POST',
                    body: JSON.stringify({
                      type: item.spec,
                      inputs,
                    })
                  })
                  const res = z.string().parse(await req.json())
                  router.push(`${session_id ? `/session/${session_id}` : ''}/graph/${res}`)
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
          </div>
        )
      }}</Catalog>
    </>
  )
}
