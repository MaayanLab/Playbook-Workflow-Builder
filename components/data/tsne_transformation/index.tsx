import python from '@/utils/python'
import { BokehPlot } from '@/components/viz/bokeh'
import { MetaNode } from '@/spec/metanode'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { scatterplot_icon } from '@/icons'

export const TSNEBokehPlotFromGeneCountMatrix = MetaNode('TSNEBokehPlotFromGeneCountMatrix')
  .meta({
    label: 'TSNE Bokeh Plot From Gene Count Matrix',
    description: 'Construct t-SNE bokeh plot From gene count matrix',
    icon: [scatterplot_icon],
  })
  .inputs({ matrix: GeneCountMatrix })
  .output(BokehPlot)
  .resolve(async (props) => await python(
      'components.data.tsne_transformation.tsne_transformation',
    { kargs: [props.inputs.matrix] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The gene count matrix was then visualized as a t-SNE plot\\ref{van der Maaten, L., & Hinton, G. (2008). Visualizing Data using t-SNE. Journal of Machine Learning Research, 9(86), 2579-2605. Retrieved from http://jmlr.org/papers/v9/vandermaaten08a.html} ${props.output_ref}.`,
    introduction: `t-Distributed Stochastic Neighbor Embedding (t-SNE) is a technique that improves on Stochastic Neighbor Embedding to produce optimized visualizations of high-dimensional data.`,
    methods: `Data in the gene count matrix are visualized using t-SNE and plotted using the Bokeh toolkit\\ref{Bokeh Development Team (2018). Bokeh: Python library for interactive visualization, http://www.bokeh.pydata.org}.`,
    legend: `A Bokeh plot displaying the t-SNE visualization of the gene count matrix.`,
  }))
  .build()
