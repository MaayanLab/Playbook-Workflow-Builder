import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneSignature } from '@/components/data/gene_signature'
import { norm_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'



export const volcanoplot = MetaNode('volcanoplot')
  .meta({
    label: 'volcanoplot',
    description: 'volcanoplot',
    icon: [norm_icon],
  })
  .inputs({ sig: GeneSignature })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.volcanoplot.createvolcano',
    { kargs: [props.inputs.sig]  },
  ))
  .story(props =>
    `The gene count matrix was then analyzed by volcano plot.`
  )
  .build()


  