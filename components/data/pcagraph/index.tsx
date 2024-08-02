import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { norm_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'
import { AnnData } from '@/components/data/anndata'


export const PCAGraphWithoutMetadata = MetaNode('PCAGraphWithoutMetadata')
  .meta({
    label: 'PCA Graph of RNA-seq Data from Gene Count Matrix',
    description: 'Construct an interactive, three-dimensional scatter plot of the first three Principal Components (PCs) of data from a Gene Count Matrix. Each point represents an RNA-seq sample. Samples with similar gene expression profiles are closer in the three-dimensional space.',
    icon: [norm_icon],
  })
  .inputs({ matrix: GeneCountMatrix })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.pcagraph.createpcanometa',
    { kargs: [props.inputs.matrix] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The gene count matrix was then visualized as a PCA plot.`,
    introduction: `Principal Component Analysis (PCA) is a linear transformation technique widely used in gene expression analysis and metabolomics studies. PCA plays a predominant role in dimensionality reduction as well as visualisation for high dimensional metabolomics data by retaining only the principal components that contain the most amount of information thereby simplifying the dataset and making it easier for downstream analytical processes.`,
    methods: `Dimensionality reduction of the gene count matrix was performed using PCA with the normalization set to log-counts-per-million (logCPM). The first three principle components (PCs) were used to generate an interactive scatter plot.`,
    legend: `An interactive scatterplot showing the first three principal components of the gene count matrix.`,
  }))
  .build()

  export const PCAGraphWithMetadata = MetaNode('PCAGraphWithMetadata')
  .meta({
    label: 'PCA Graph of RNA-seq Data from AnnData file',
    description: 'Construct an interactive, three-dimensional scatter plot of the first three Principal Components (PCs) of data from an AnnData file. Each point represents an RNA-seq sample. Samples with similar gene expression profiles are closer in the three-dimensional space.',
    icon: [norm_icon],
  })
  .inputs({ anndata: AnnData })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.pcagraph.createmetapcagraph',
    { kargs: [props.inputs.anndata] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The AnnData file was then visualized as a PCA graph.`,
    introduction: `Principal Component Analysis (PCA) is a linear transformation technique widely used in gene expression analysis and metabolomics studies. PCA plays a predominant role in dimensionality reduction as well as visualisation for high dimensional metabolomics data by retaining only the principal components that contain the most amount of information thereby simplifying the dataset and making it easier for downstream analytical processes.`,
    methods: `Metadata was uploaded from an AnnData file. The AnnData Python package handles annotated data matrices in memory and on disk\\ref{doi:10.1101/2021.12.16.473007}. Dimensionality reduction of the metadata was performed using PCA with the normalization set to log-counts-per-million (logCPM). The first three principle components (PCs) were used to generate an interactive scatter plot.`,
    legend: `An interactive scatterplot showing the first three principal components of the AnnData.`,
  }))
  .build()