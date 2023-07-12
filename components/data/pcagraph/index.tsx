import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { norm_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'



export const PCAgraph = MetaNode('PCAgraph')
  .meta({
    label: 'PCA graph of RNAseq data',
    description: 'analyze gene count matrix of RNAseq data and return PCA graph',
    icon: [norm_icon],
  })
  .inputs({ matrix: GeneCountMatrix })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.pcagraph.createpca',
    { kargs: [props.inputs.matrix]  },
  ))
  .story(props =>
    `The gene count matrix was then analyzed by PCA. The PCA plug-in embeds an interactive, three-dimensional scatter plot of the first three Principal Components (PCs) of the data. Each point represents an RNA-seq sample. Samples with similar gene expression profiles are closer in the three-dimensional space.`
  )
  .build()