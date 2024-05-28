import python from '@/utils/python'
import { PlotlyPlot } from '@/components/viz/plotly'
import { MetaNode } from '@/spec/metanode'
import { plot_icon } from '@/icons'
import { GMT } from '@/components/data/gene_matrix_transpose'
import { DMT } from '@/components/data/drug_matrix_transpose'

export const UMAPFromXMT = [GMT, DMT].map(XMT =>
  MetaNode(`UMAPFrom[${XMT.spec}]`)
    .meta({
      label: `UMAP plot from ${XMT.meta.label}`,
      description: `Construct UMAP plot with ${XMT.meta.label}`,
      icon: [plot_icon],
    })
    .inputs({ matrix: XMT })
    .output(PlotlyPlot)
    .resolve(async (props) => await python(
      'components.viz.idf_umap.idf_umap',
      { kargs: [props.inputs.matrix] },
      message => props.notify({ type: 'info', message }),
    ))
    .story(props =>
      `UMAP was applied to the inverse document frequency of the ${XMT.meta.label.toLocaleLowerCase()}.`
    )
    .build()
)
