import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { Gene, Drug } from '@/components/core/input/primitives'

const Set_T = (T: typeof Gene) => MetaNode.createData(`Set[${T.name}]`)
  .meta({
    label: `${T.label} Set`,
    description: `Set of ${T.label}s`,
    color: T.color,
    example: T.examples.set,
  })
  .codec(z.array(z.string()))
  .view(set => {
    return <textarea readOnly value={set.join('\n')} />
  })
  .build()

export const GeneSet = Set_T(Gene)
export const DrugSet = Set_T(Drug)

const Input_Set_T = (T: typeof GeneSet) => MetaNode.createProcess(`Input[${T.spec}]`)
  .meta({
    label: `${T.meta.label} Set Input`,
    description: `Start with a set of ${T.meta.label}`,
  })
  .inputs()
  .output(T)
  .prompt(props => {
    const [set, setSet] = React.useState('')
    React.useEffect(() => { setSet((props.output||[]).join('\n')) }, [props.output])
    return (
      <div>
        <textarea value={set} onChange={evt => setSet(evt.target.value)} />
        <button onClick={evt => props.submit(set.split(/\r?\n/g))}>Submit</button>
        <button onClick={evt => props.submit(T.meta.example)}>Example</button>
      </div>
    )
  })
  .build()

export const InputGeneSet = Input_Set_T(GeneSet)
export const InputDrugSet = Input_Set_T(DrugSet)
