import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneSignature } from '@/components/data/gene_signature'
import { norm_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'



export const differentialexpression = MetaNode('differentialexpression')
  .meta({
    label: 'differentialexpression',
    description: 'differentialexpression',
    icon: [norm_icon],
  })
  .inputs({ Matrix: GeneSignature })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.differentialexpression.createdifferentialexpression',
    { kargs: [props.inputs.Matrix]  },
  ))
  .story(props =>
    `differentialexpression.`
  )
  .build()


