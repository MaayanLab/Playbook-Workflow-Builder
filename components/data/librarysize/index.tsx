import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { norm_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'



export const VisualizeLibrarySizes = MetaNode('VisualizeLibrarySizes')
  .meta({
    label: 'Library Size Bar Plot from Gene Count Matrix',
    description: 'Construct a bar plot which displays the total number of reads mapped to each RNA-seq sample in the dataset',
    icon: [norm_icon],
  })
  .inputs({ matrix: GeneCountMatrix })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.librarysize.createlibrarysize',
    { kargs: [props.inputs.matrix]  },
  ))
  .story(props =>
    `The gene count matrix was then visualized as a bar plot representing library sizes.`
  )
  .build()


