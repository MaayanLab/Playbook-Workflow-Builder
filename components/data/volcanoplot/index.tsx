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
    `The gene count matrix was then analyzed by volcano plot. The volcano plot plug-in embeds a scatter plot which displays the log2-fold changes and statistical significance of each gene calculated by performing differential gene expression analysis comparing samples in the Control group to samples in the Perturbation group. Every point in the plot represents a gene.`
  )
  .build()


  