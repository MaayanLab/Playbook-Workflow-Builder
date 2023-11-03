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
  .story(props =>
    `The gene count matrix was then visualized as a t-SNE plot${''/* [FIGURE]*/}.`
  )
  .build()
