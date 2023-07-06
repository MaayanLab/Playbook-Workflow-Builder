import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { norm_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'



export const librarysize = MetaNode('librarysize')
  .meta({
    label: 'librarysize',
    description: 'librarysize',
    icon: [norm_icon],
  })
  .inputs({ matrix: GeneCountMatrix })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.librarysize.createlibrarysize',
    { kargs: [props.inputs.matrix]  },
  ))
  .story(props =>
    `librarysize.`
  )
  .build()


