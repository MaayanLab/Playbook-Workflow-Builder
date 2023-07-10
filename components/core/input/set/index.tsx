import React from 'react'
import { DataMetaNode, InternalDataMetaNode, MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import {
  Disease,
  Drug,
  Gene,
  Metabolite,
  Pathway,
  Phenotype,
  Primative,
  RegulatoryElement,
  Tissue,
  Variant,
} from '@/components/core/primitives'
import { Table, Cell, Column } from '@/app/components/Table'
import { input_icon, set_icon } from '@/icons'
import * as array from '@/utils/array'
import { downloadBlob } from '@/utils/download'
import dynamic from 'next/dynamic'
import pluralize from 'pluralize'

const Bp4Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))
const Bp4TextArea = dynamic(() => import('@blueprintjs/core').then(({ TextArea }) => TextArea))

const Set_T = (T: Primative) => MetaNode(`Set[${T.name}]`)
  .meta({
    label: `${T.label} Set`,
    description: `Set of ${pluralize(T.label)}`,
    icon: [...array.ensureArray(T.icon), set_icon],
    color: T.color,
    tags: {
      Type: {
        [T.label]: 1,
      },
      Cardinality: {
        Set: 1,
      },
    },
    ...(T.extra?.set?.meta || {}),
  })
  .codec(z.object({ description: z.string().optional(), set: z.array(z.string()) }))
  .view(({ set }) => {
    return (
      <Table
        height={500}
        cellRendererDependencies={[set]}
        numRows={set.length}
        enableGhostCells
        enableFocusedCell
        downloads={{
          JSON: () => downloadBlob(new Blob([JSON.stringify(set)], { type: 'application/json;charset=utf-8' }), 'data.json'),
          GMT: () => downloadBlob(new Blob([[`${T.label} Set`, '', ...set].join('\t')], { type: 'text/tab-separated-values;charset=utf-8' }), 'data.gmt'),
        }}
      >
        <Column
          name={T.label}
          cellRenderer={row => <Cell key={row+''}>{set[row]}</Cell>}
        />
      </Table>
    )
  })
  .build()

export const DiseaseSet = Set_T(Disease)
export const DrugSet = Set_T(Drug)
export const GeneSet = Set_T(Gene)
export const VariantSet = Set_T(Variant)
export const RegulatoryElementSet = Set_T(RegulatoryElement)
export const PathwaySet = Set_T(Pathway)
export const PhenotypeSet = Set_T(Phenotype)
export const TissueSet = Set_T(Tissue)
export const MetaboliteSet = Set_T(Metabolite)

const Input_Set_T = (T: Primative, SetT: DataMetaNode<InternalDataMetaNode & { data: { description?: string, set: string[] } }>) => MetaNode(`Input[${SetT.spec}]`)
  .meta({
    label: `${T.label} Set Input`,
    description: `Start with a set of ${pluralize(T.label)}`,
    icon: [input_icon],
    tags: {
      Type: {
        [T.label]: 1,
      },
      Cardinality: {
        Set: 1,
      },
    },
    ...(T.extra?.set?.meta || {}),
  })
  .inputs()
  .output(SetT)
  .prompt(props => {
    const [set, setSet] = React.useState('')
    const [description, setDescription] = React.useState('')
    React.useEffect(() => {
      setSet(((props.output||{}).set||[]).join('\n'))
      setDescription(((props.output||{}).description||''))
    }, [props.output])
    return (
      <div>
        <Bp4TextArea
          placeholder="Newline separated set of terms"
          rows={8}
          fill
          large
          onChange={evt => setSet(evt.target.value)}
          value={set}
        />
        <div className="bp4-input-group">
          <input
            type="text"
            className="bp4-input"
            placeholder={`${T.label} Set description`}
            onChange={evt => setDescription(evt.target.value)}
            value={description}
          />
        </div>
        {T.extra?.set?.meta?.example !== undefined ?
          <Bp4Button
            large
            rightIcon="send-to-graph"
            onClick={evt => {
              if (T.extra?.set?.meta?.example !== undefined) {
                setDescription(T.extra.set.meta.example.description)
                setSet(T.extra.set.meta.example.set.join('\n'))
              }
            }}
            text="Example"
          />
          : null}
        <Bp4Button
          large
          type="submit"
          text="Submit"
          rightIcon="bring-data"
          onClick={evt => props.submit({ description: description, set: set.split(/\r?\n/g).filter(line => line.replace(/^\s+/, '').replace(/\s+$/, '')) })}
          disabled={set.length === 0}
        />
      </div>
    )
  })
  .story(props =>
    `The workflow starts with a ${T.label.toLocaleLowerCase()} set${props.output && props.output.description ? ` created from ${props.output.description}` : ''}.`
  )
  .build()

export const InputGeneSet = Input_Set_T(Gene, GeneSet)
export const InputVariantSet = Input_Set_T(Variant, VariantSet)
export const InputRegulatoryElementSet = Input_Set_T(RegulatoryElement, RegulatoryElementSet)
export const InputDrugSet = Input_Set_T(Drug, DrugSet)
export const InputMetaboliteSet = Input_Set_T(Metabolite, MetaboliteSet)

import {
  DrugTerm,
  GeneTerm,
  MetaboliteTerm,
  VariantTerm,
} from '@/components/core/input/term'
export const TermToSetT = [
  { SetT: GeneSet, TermT: GeneTerm },
  { SetT: VariantSet, TermT: VariantTerm },
  { SetT: DrugSet, TermT: DrugTerm },
  { SetT: MetaboliteSet, TermT: MetaboliteTerm },
].map(({ SetT, TermT }) =>
  MetaNode(`TermToSet[${TermT.spec}]`)
    .meta({
      label: `Create Set from ${pluralize(TermT.meta.label)}`,
      description: 'Join several terms into a set'
    })
    .inputs({ terms: [TermT] })
    .output(SetT)
    .resolve(async (props) => {
      return { set: props.inputs.terms }
    })
    .story(props =>
      `The specified genes were combined into one gene set.`
    )
    .build()
)
