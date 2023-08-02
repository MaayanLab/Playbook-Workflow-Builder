import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { norm_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'
import { AnnData } from '@/components/data/anndata'



export const VisualizeLibrarySizes = MetaNode('VisualizeLibrarySizes')
  .meta({
    label: 'Library Size Bar Plot from Gene Count Matrix',
    description: 'Construct a bar plot which displays the total number of reads mapped to each RNA-seq sample in the dataset from a gene count matrix.',
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


  export const VisualizeLibrarySizesfromAnnData = MetaNode('VisualizeLibrarySizesFromAnnData')
  .meta({
    label: 'Library Size Bar Plot from AnnData File',
    description: 'Construct a bar plot which displays the total number of reads mapped to each RNA-seq sample in the dataset from an AnnData file.',
    icon: [norm_icon],
  })
  .inputs({
    anndata: AnnData,
  })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.librarysize.createlibrarysizefromanndata',
    { kargs: [props.inputs.anndata] },
  ))
  .story(props =>
    `The AnnData file was then visualized as a bar plot representing library sizes.`
  )
  .build()



