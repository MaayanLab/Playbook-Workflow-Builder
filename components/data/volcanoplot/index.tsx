import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneSignature } from '@/components/data/gene_signature'
import { norm_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'



export const VolcanoPlot = MetaNode('VolcanoPlot')
  .meta({
    label: 'Volcano Plot from Differential Expression Table',
    description: 'Construct a scatter plot which displays the log2-fold changes and statistical significance of each gene calculated by performing differential gene expression analysis comparing samples in the Control group to samples in the Perturbation group. Every point in the plot represents a gene.',
    icon: [norm_icon],
  })
  .inputs({ sig: GeneSignature })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.volcanoplot.createvolcano',
    { kargs: [props.inputs.sig]  },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The data in the differential expression table was then visualized as a volcano plot.`,
  }))
  .build()


  