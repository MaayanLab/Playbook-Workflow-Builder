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
  ))
  .story(props =>
    `The gene count matrix was then visualized as a UMAP plot${''/* [FIGURE]*/} [\\ref{doi:10.48550/arXiv.1802.03426}].`
  )
  .build()
