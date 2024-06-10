import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { MetaboliteCountMatrix } from '@/components/data/metabolite_count_matrix'
import { norm_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'
import { MetAnnData } from '@/components/data/metanndata'


export const MetPCAGraphWithNoMeta = MetaNode('MetPCAGraphWithoutMetadata')
  .meta({
    label: 'PCA Graph of Metabolomics Data from Metabolite Count Matrix',
    description: 'Construct an interactive, three-dimensional scatter plot of the first three Principal Components (PCs) of data from a Metabolite Count Matrix. Each point represents an Metabolomics sample. Samples with similar metabolite expression profiles are closer in the three-dimensional space.',
    icon: [norm_icon],
  })
  .inputs({ matrix: MetaboliteCountMatrix })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.metpcagraph.createmetpcanometa',
    { kargs: [props.inputs.matrix] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The metabolite count matrix was then visualized as a PCA plot.`,
  }))
  .build()

export const MetPCAGraphWithMetadata = MetaNode('MetPCAGraphWithMetadata')
  .meta({
    label: 'PCA Graph of Metabolite Data from MetAnnData file',
    description: 'Construct an interactive, three-dimensional scatter plot of the first three Principal Components (PCs) of data from an MetAnnData file. Each point represents an metabolomics sample. Samples with similar metabolite expression profiles are closer in the three-dimensional space.',
    icon: [norm_icon],
  })
  .inputs({ anndata: MetAnnData })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.metpcagraph.createmetametpcagraph',
    { kargs: [props.inputs.anndata] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The MetAnnData file for metabolites was then visualized as a PCA graph.`,
  }))
  .build()