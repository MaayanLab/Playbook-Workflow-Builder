import React from 'react'
import { MetaNode } from '@/spec/metanode'
import * as t from 'io-ts'
import codecFrom from '@/utils/io-ts-codec'

const example = 'ACE2'

export const GeneSymbol = MetaNode.createData('GeneSymbol')
  .meta({
    label: 'GeneSymbol',
    description: 'An unresolved Gene Symbol',
    example,
  })
  .codec(codecFrom(t.string))
  .view(gene => (
    <div>{gene} (gene)</div>
  ))
  .build()

export const GeneSymbolInput = MetaNode.createProcess('GeneSymbolInput')
  .meta({
    label: 'Input a Gene',
    description: 'A gene input prompt',
    default: '',
  })
  .inputs()
  .output(GeneSymbol)
  .prompt(props => {
    const [gene, setGene] = React.useState('')
    React.useEffect(() => { setGene(props.output) }, [props.output])
    return (
      <div>
        <input value={gene} onChange={evt => setGene(evt.target.value)} />
        <button onClick={evt => props.submit(gene)}>Submit</button>
        <button onClick={evt => props.submit(example)}>Example</button>
      </div>
    )
  })
  .build()
