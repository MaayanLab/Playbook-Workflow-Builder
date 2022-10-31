import python from '@/utils/python'
import { PlotlyPlot } from '@/components/viz/plotly'
import { MetaNode } from '@/spec/metanode'
import { SignificantTissues } from '@/components/core/significant_tissues'
import { barchart_icon } from '@/icons'

export const TissueBarplotFromSignificantTissue = MetaNode.createProcess('TissueBarplotFromSignificantTissue')
  .meta({
    label: 'Tissue Barplot from Significant Tissues',
    description: 'Construct Tissue Barplot with Significant Tissues',
    icon: [barchart_icon],
  })
  .codec()
  .inputs({ tissues: SignificantTissues })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.viz.tissue_barplot.tissue_barplot',
    { kargs: [props.inputs.tissues] },
  ))
  .build()
