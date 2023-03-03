import python from '@/utils/python'
import { PlotlyPlot } from '@/components/viz/plotly'
import { MetaNode } from '@/spec/metanode'
import { ScoredTissues } from '@/components/core/input/scored'
import { barchart_icon } from '@/icons'

export const TissueBarplotFromSignificantTissue = MetaNode('TissueBarplotFromSignificantTissue')
  .meta({
    label: 'Tissue Barplot from Significant Tissues',
    description: 'Construct Tissue Barplot with Significant Tissues',
    icon: [barchart_icon],
  })
  .inputs({ tissues: ScoredTissues })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.viz.tissue_barplot.tissue_barplot',
    { kargs: [props.inputs.tissues] },
  ))
  .story(props =>
    `To visualize the level of expression across tissues, a bar plot was created [FIGURE].`
  )
  .build()
