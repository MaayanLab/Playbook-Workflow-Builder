import python from '@/utils/python'
import { BokehPlot } from '@/components/viz/bokeh'
import { MetaNode } from '@/spec/metanode'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { scatterplot_icon } from '@/icons'

export const UMAPBokehPlotFromGeneCountMatrix = MetaNode('UMAPBokehPlotFromGeneCountMatrix')
  .meta({
    label: 'UMAP Bokeh Plot From Gene Count Matrix',
    description: 'Construct UMAP bokeh plot From gene count matrix',
    icon: [scatterplot_icon],
  })
  .inputs({ matrix: GeneCountMatrix })
  .output(BokehPlot)
  .resolve(async (props) => await python(
      'components.data.umap_transformation.umap_transformation',
    { kargs: [props.inputs.matrix]  },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The gene count matrix was then visualized as a UMAP plot\\ref{doi:10.48550/arXiv.1802.03426} ${props.output_ref}.`,
    introduction: `Uniform Manifold Approximation and Projection (UMAP) is a dimensionality reduction technique that is useful to visualize high-dimensional data\\ref{doi:10.48550/arXiv.1802.03426}.`,
    methods: `A UMAP representation was generated for the gene count matrix and plotted using the Bokeh toolkit\\ref{Bokeh Development Team (2018). Bokeh: Python library for interactive visualization, http://www.bokeh.pydata.org}.`,
    legend: `A Bokeh plot displaying the UMAP representation of the gene count matrix.`,
  }))
  .build()
