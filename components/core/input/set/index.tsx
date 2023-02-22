import React from 'react'
import { DataMetaNode, InternalDataMetaNode, MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { Gene, Drug, Primative, Pathway, Phenotype, Tissue, Disease } from '@/components/core/input/primitives'
import { Table, Cell, Column } from '@/app/components/Table'
import { input_icon, set_icon } from '@/icons'
import * as array from '@/utils/array'
import { downloadBlob } from '@/utils/download'
import dynamic from 'next/dynamic'

const Bp4Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))
const Bp4TextArea = dynamic(() => import('@blueprintjs/core').then(({ TextArea }) => TextArea))

const Set_T = (T: Primative) => MetaNode(`Set[${T.name}]`)
  .meta({
    label: `${T.label} Set`,
    description: `Set of ${T.label}s`,
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
  .codec(z.array(z.string()))
  .view(set => {
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
export const PathwaySet = Set_T(Pathway)
export const PhenotypeSet = Set_T(Phenotype)
export const TissueSet = Set_T(Tissue)

const Input_Set_T = (T: Primative, SetT: DataMetaNode<InternalDataMetaNode & { data: string[] }>) => MetaNode(`Input[${SetT.spec}]`)
  .meta({
    label: `${T.label} Set Input`,
    description: `Start with a set of ${T.label}s`,
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
    React.useEffect(() => { setSet((props.output||[]).join('\n')) }, [props.output])
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
        {T.extra?.set?.meta?.example !== undefined ?
          <Bp4Button
            large
            rightIcon="send-to-graph"
            onClick={evt => {
              if (T.extra?.set?.meta?.example !== undefined) {
                setSet(T.extra.set.meta.example.join('\n'))
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
          onClick={evt => props.submit(set.split(/\r?\n/g))}
        />
      </div>
    )
  })
  .build()

export const InputGeneSet = Input_Set_T(Gene, GeneSet)
export const InputDrugSet = Input_Set_T(Drug, DrugSet)
