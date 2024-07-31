import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { MetaboliteCountMatrix } from '@/components/data/metabolite_count_matrix'
import { norm_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'
import { MetAnnData } from '@/components/data/metanndata'


export const MetPCAGraphWithoutMetadata = MetaNode('MetPCAGraphWithoutMetadata')
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
    introduction: `Principal Component Analysis is a linear transformation technique widely used in gene expression analysis and metabolomics studies. PCA plays a predominant role in dimensionality reduction as well as visualisation for high dimensional metabolomics data by retaining only the principal components that contain the most amount of information thereby simplifying the dataset and making it easier for downstream analytical processes.`,
    methods: `A log counts per million (logCPM) transformation and Principal Component Analysis (PCA) is performed on metabolomics count matrix. The \`logcpm_transform\` function, implemented with a \`FunctionTransformer\` using \`np.log1p\`, applies a log transformation to the dataset. This step normalizes the metabolite abundance data and reduces skewness, making the data more suitable for analysis. The transformed data, initially in a NumPy array, is then converted into a pandas DataFrame for easier manipulation and further analysis. Finally, PCA is conducted to reduce the dimensionality of the transformed metabolomics data to three principal components, which helps in identifying patterns and visualizing the data's underlying structure.`,
    legend: `The results are visualized as a 3D scatter plot (with the ability to pan and rotate along the three axes) with the three prinicpial components along the three axes.`,
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
    introduction: `Principal Component Analysis is a linear transformation technique widely used in gene expression analysis and metabolomics studies. PCA plays a predominant role in dimensionality reduction as well as visualisation for high dimensional metabolomics data by retaining only the principal components that contain the most amount of information thereby simplifying the dataset and making it easier for downstream analytical processes.`,
    methods: `The presence of annotations helps clearly differentiate between the control and perturbation groups in the metabolite count matrix, providing a visual legend for the scatter plot. Annotated metabolomics data is read, including metadata that distinguishes between control and perturbation groups. Relevant metadata tis then extracted to identify control and perturbed case samples, creating masks to distinguish them. PCA is run on the dataset, transposed for correct dimensional alignment and a 3D scatter plot using Plotly, with points colored blue for control samples and red for perturbation samples. Axis labels showing the variance explained by each principal component, with appropriate font sizing for readability.`,
    legend: `The results are visualized as an annotated 3D scatter plot (with the ability to pan and rotate along the three axes) with the three prinicpial components along the three axes with control and perturbed case distinguished and annotated.`,
  }))
  .build()