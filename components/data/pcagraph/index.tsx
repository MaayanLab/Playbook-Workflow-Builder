import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { norm_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'



export const PCAGraph = MetaNode('PCAGraph')
  .meta({
    label: 'PCA Graph of RNAseq Data from Gene Count Matrix',
    description: 'Construct an interactive, three-dimensional scatter plot of the first three Principal Components (PCs) of data from a Gene Count Matrix. Each point represents an RNA-seq sample. Samples with similar gene expression profiles are closer in the three-dimensional space.',
    icon: [norm_icon],
  })
  .inputs({ matrix: GeneCountMatrix })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.pcagraph.createpca',
    { kargs: [props.inputs.matrix]  },
  ))
  .story(props =>
    `The gene count matrix was then visualized as a PCA graph.`
  )
  .build()

  