import React from 'react'
import { MetaNode, MetaNodeDataType } from '@/spec/metanode'
import { z } from 'zod'
import { Intent } from '@blueprintjs/core'
import dynamic from 'next/dynamic'

const TextArea = dynamic(() => import('@blueprintjs/core').then(({ TextArea }) => TextArea))

export const Suggestion = MetaNode.createData('Suggestion')
  .meta({
    label: 'Suggestion',
    description: 'An actual suggestion',
  })
  .codec(z.string())
  .view(suggestion => <div>{suggestion}</div>)
  .build()

export const SuggestInputEdge = MetaNode.createProcess('SuggestInputEdge')
  .meta({
    label: 'Suggest an input type',
    description: 'Provide a description about what should be here',
  })
  .inputs()
  .output(Suggestion)
  .prompt((props) => {
    const [suggestion, setSuggestion] = React.useState('')
    return (
      <TextArea
        growVertically={true}
        large={true}
        intent={Intent.PRIMARY}
        onChange={(evt) => setSuggestion(evt.target.value)}
        value={suggestion}
      />
    )
  })
  .build()


import { FileURL } from '@/components/core/file'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { GeneTerm } from '@/components/core/input/term'
import { GeneInfo } from '@/components/service/mygeneinfo'
import { SignificantTissues } from '@/components/core/significant_tissues'

export const SuggestInteractiveEdge = [
  FileURL,
  GeneCountMatrix,
  GeneTerm,
  GeneInfo,
  SignificantTissues,
].flatMap(input => [
  MetaNode.createProcess(`SuggestInteractiveEdge[${input.spec}]`)
    .meta({
      label: 'Suggest a visualization method',
      description: `This would visualize the ${input.meta.label || input.spec}. Provide a description about what should be here.`,
    })
    .inputs({ input } as Record<string, MetaNodeDataType>)
    .output(Suggestion)
    .prompt((props) => {
      const [suggestion, setSuggestion] = React.useState('')
      return (
        <TextArea
          growVertically={true}
          large={true}
          intent={Intent.PRIMARY}
          onChange={(evt) => setSuggestion(evt.target.value)}
          value={suggestion}
        />
      )
    })
    .build(),
  MetaNode.createProcess(`SuggestResolveEdge[${input.spec}]`)
    .meta({
      label: 'Suggest an algorithm or data transformation method',
      description: `This would transform the ${input.meta.label || input.spec}. Provide a description about what should be here.`,
    })
    .inputs({ input } as Record<string, MetaNodeDataType>)
    .output(Suggestion)
    .prompt((props) => {
      const [suggestion, setSuggestion] = React.useState('')
      return (
        <TextArea
          growVertically={true}
          large={true}
          intent={Intent.PRIMARY}
          onChange={(evt) => setSuggestion(evt.target.value)}
          value={suggestion}
        />
      )
    })
    .build()
  ]
)
