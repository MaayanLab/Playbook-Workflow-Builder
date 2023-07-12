import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { norm_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'



export const VisualizeLibrarySize = MetaNode('VisualizeLibrarySize')
  .meta({
    label: 'VisualizeLibrarySize',
    description: 'VisualizeLibrarySize',
    icon: [norm_icon],
  })
  .inputs({ matrix: GeneCountMatrix })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.librarysize.createlibrarysize',
    { kargs: [props.inputs.matrix]  },
  ))
  .story(props =>
    `This library size plug-in embeds an interactive bar chart which displays the total number of reads mapped to each RNA-seq sample in the dataset.`
  )
  .build()


