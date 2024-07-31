import python from '@/utils/python'
import { BokehPlot } from '@/components/viz/bokeh'
import { MetaNode } from '@/spec/metanode'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { scatterplot_icon } from '@/icons'

export const PCABokehPlotFromGeneCountMatrix = MetaNode('PCABokehPlotFromGeneCountMatrix')
  .meta({
    label: 'PCA Bokeh Plot From Gene Count Matrix',
    description: 'Construct PCA Bokeh Plot From Gene Count Matrix',
    icon: [scatterplot_icon],
  })
  .inputs({ matrix: GeneCountMatrix })
  .output(BokehPlot)
  .resolve(async (props) => await python(
      'components.data.pca_transformation.pca_transformation',
    { kargs: [props.inputs.matrix] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The gene count matrix was then visualized as a PCA plot ${''/*({props.output_ref})*/}.`,
    introduction: `Principal Component Analysis (PCA) is a linear transformation technique widely used in gene expression analysis and metabolomics studies. PCA plays a predominant role in dimensionality reduction as well as visualisation for high dimensional metabolomics data by retaining only the principal components that contain the most amount of information thereby simplifying the dataset and making it easier for downstream analytical processes.`,
    methods: `Dimensionality reduction of the gene count matrix was performed using PCA with the normalization set to log-counts-per-million (logCPM). The results were plotted using the Bokeh toolkit\\ref{https://bokeh.org}.`,
    legend: `A Bokeh plot of PCA results for the gene count matrix.`,
  }))
  .build()
