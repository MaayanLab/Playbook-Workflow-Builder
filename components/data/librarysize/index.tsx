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
    { kargs: [props.inputs.matrix] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The gene count matrix was then visualized as a bar plot representing library sizes.`,
    introduction: `The first step in determining gene expression levels is to align reads to the sample genome. Tracking the number of reads mapped to each sample and the genes they represent is a method of quantifying gene expression.`,
    methods: `Expression data was quantified as gene-level counts by mapping transcript sequences to genes using the ARCHS4 pipeline\\ref{doi:10.1038/s41467-018-03751-6}.`,
    legend: `A Plotly bar plot representing library sizes.`,
  }))
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
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The AnnData file was then visualized as a bar plot representing library sizes.`,
    introduction: `The first step in determining gene expression levels is to align reads to the sample genome. Tracking the number of reads mapped to each sample and the genes they represent is a method of quantifying gene expression.  The AnnData Python package handles annotated data matrices in memory and on disk\\ref{doi:10.1101/2021.12.16.473007} and can be used to store gene count information.`,
    methods: `Expression data was quantified as gene-level counts by mapping transcript sequences to genes using the ARCHS4 pipeline\\ref{doi:10.1038/s41467-018-03751-6}.`,
    legend: `A Plotly bar plot representing library sizes.`,
  }))
  .build()



