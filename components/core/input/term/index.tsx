import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { Gene, Drug } from '@/components/core/input/primitives'

const Term_T = (T: typeof Gene) => MetaNode.createData(`Term[${T.name}]`)
  .meta({
    label: T.label,
    description: `${T.label} Term`,
    color: T.color,
    example: T.examples.term,
  })
  .codec(z.string())
  .view(term => {
    return <div>{T.label}: {term}</div>
  })
  .build()

export const GeneTerm = Term_T(Gene)
export const DrugTerm = Term_T(Drug)

const Input_Term_T = (T: typeof GeneTerm) => MetaNode.createProcess(`Input[${T.spec}]`)
  .meta({
    label: `${T.meta.label} Input`,
    description: `Start with a ${T.meta.label} term`,
  })
  .inputs()
  .output(T)
  .prompt(props => {
    const [term, setTerm] = React.useState('')
    React.useEffect(() => { setTerm(props.output||'') }, [props.output])
    return (
      <div>
        <input value={term} onChange={evt => setTerm(evt.target.value)} />
        <button onClick={evt => props.submit(term)}>Submit</button>
        <button onClick={evt => props.submit(T.meta.example)}>Example</button>
      </div>
    )
  })
  .build()

export const InputGeneTerm = Input_Term_T(GeneTerm)
export const InputDrugTerm = Input_Term_T(DrugTerm)
