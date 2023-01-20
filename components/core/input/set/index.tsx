import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { Gene, Drug, Primative, Pathway, Phenotype, Tissue, Disease } from '@/components/core/input/primitives'
import { Button, TextArea } from '@blueprintjs/core'
import { Table2 as Table, Column, Cell } from '@blueprintjs/table'
import { input_icon, set_icon } from '@/icons'

const Set_T = (T: Primative) => MetaNode.createData(`Set[${T.name}]`)
  .meta({
    label: `${T.label} Set`,
    description: `Set of ${T.label}s`,
    icon: [...T.icon, set_icon],
    color: T.color,
    example: T.examples.set,
  })
  .codec(z.array(z.string()))
  .view(set => {
    return (
      <div style={{ height: 500 }}>
        <Table
          cellRendererDependencies={[set]}
          numRows={set.length}
          enableGhostCells
          enableFocusedCell
        >
          <Column
            name={T.label}
            cellRenderer={row => <Cell key={row+''}>{set[row]}</Cell>}
          />
        </Table>
      </div>
    )
  })
  .build()

export const DiseaseSet = Set_T(Disease)
export const DrugSet = Set_T(Drug)
export const GeneSet = Set_T(Gene)
export const PathwaySet = Set_T(Pathway)
export const PhenotypeSet = Set_T(Phenotype)
export const TissueSet = Set_T(Tissue)

const Input_Set_T = (T: typeof GeneSet) => MetaNode.createProcess(`Input[${T.spec}]`)
  .meta({
    label: `${T.meta.label} Input`,
    description: `Start with a set of ${T.meta.label}`,
    icon: [input_icon],
  })
  .inputs()
  .output(T)
  .prompt(props => {
    const [set, setSet] = React.useState('')
    React.useEffect(() => { setSet((props.output||[]).join('\n')) }, [props.output])
    return (
      <div>
        <TextArea
          placeholder="Newline separated set of terms"
          rows={8}
          fill
          large
          onChange={evt => setSet(evt.target.value)}
          value={set}
        />
        <Button large rightIcon="bring-data" onClick={evt => props.submit(set.split(/\r?\n/g))}>Submit</Button>
        {T.meta.example !== undefined ?
          <Button
            large
            rightIcon="send-to-graph"
            onClick={evt => {
              if (T.meta.example !== undefined) {
                props.submit(T.meta.example)
              }
            }}>Example</Button>
          : null}
      </div>
    )
  })
  .build()

export const InputGeneSet = Input_Set_T(GeneSet)
export const InputDrugSet = Input_Set_T(DrugSet)
