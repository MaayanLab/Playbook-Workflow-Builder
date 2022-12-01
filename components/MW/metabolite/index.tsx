import React from 'react'
import { MetaNode } from '@/spec/metanode'
import * as t from 'io-ts'
import codecFrom from '@/utils/io-ts-codec'

const example = 'Glucose'

export const MetaboliteName = MetaNode.createData('MetaboliteName')
  .meta({
    label: 'Metabolite name',
    description: 'An unresolved Metabolite',
    example,
  })
  .codec(codecFrom(t.string))
  .view(metabolite => (
    <div>{metabolite} (metabolite)</div>
  ))
  .build()

export const MetaboliteNameInput = MetaNode.createProcess('MetaboliteNameInput')
  .meta({
    label: 'Input a Metabolite',
    description: 'A metabolite input prompt',
  })
  .inputs()
  .output(MetaboliteName)
  .prompt(props => {
    const [metabolite, setMetabolite] = React.useState(props.output || '')
    return (
      // https://www.metabolomicsworkbench.org/rest/refmet/name/Glucose/all/
      <div>
        <input value={metabolite} onChange={evt => setMetabolite(evt.target.value)} />
        <button onClick={evt => props.submit(metabolite)}>Submit</button>
        <button onClick={evt => props.submit(example)}>Example</button>
      </div>
    )
  })
  .build()
