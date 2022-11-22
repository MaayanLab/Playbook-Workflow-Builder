import python from '@/utils/python'
import { BokehPlot } from '@/components/bokeh'
import { MetaNode } from '@/spec/metanode'
import { GeneCountMatrix } from '@/components/gene_count_matrix'

export const TSNEBokehPlotFromGeneCountMatrix = MetaNode.createProcess('TSNEBokehPlotFromGeneCountMatrix')
  .meta({
    label: 'TSNE Bokeh Plot From Gene Count Matrix',
    description: 'Construct t-SNE bokeh plot From gene count matrix',
  })
  .codec()
  .inputs({ matrix: GeneCountMatrix })
  .output(BokehPlot)
  .resolve(async (props) => await python(
      'components.tsne_transformation.tsne_transformation',
    { kargs: [props.inputs.matrix]  },
  ))
  .build()
