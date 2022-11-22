import python from '@/utils/python'
import { BokehPlot } from '@/components/bokeh'
import { MetaNode } from '@/spec/metanode'
import { GeneCountMatrix } from '@/components/gene_count_matrix'

export const UMAPBokehPlotFromGeneCountMatrix = MetaNode.createProcess('UMAPBokehPlotFromGeneCountMatrix')
  .meta({
    label: 'UMAP Bokeh Plot From Gene Count Matrix',
    description: 'Construct UMAP bokeh plot From gene count matrix',
  })
  .codec()
  .inputs({ matrix: GeneCountMatrix })
  .output(BokehPlot)
  .resolve(async (props) => await python(
      'components.umap_transformation.umap_transformation',
    { kargs: [props.inputs.matrix]  },
  ))
  .build()
