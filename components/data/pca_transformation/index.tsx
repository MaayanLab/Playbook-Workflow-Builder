import python from '@/utils/python'
import { BokehPlot } from '@/components/viz/bokeh'
import { MetaNode } from '@/spec/metanode'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'

export const PCABokehPlotFromGeneCountMatrix = MetaNode.createProcess('PCABokehPlotFromGeneCountMatrix')
  .meta({
    label: 'PCA Bokeh Plot From Gene Count Matrix',
    description: 'Construct PCA Bokeh Plot From Gene Count Matrix',
  })
  .codec()
  .inputs({ matrix: GeneCountMatrix })
  .output(BokehPlot)
  .resolve(async (props) => await python(
      'components.data.pca_transformation.pca_transformation',
    { kargs: [props.inputs.matrix]  },
  ))
  .build()
