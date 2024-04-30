import python from '@/utils/python'
import { PlotlyPlot } from '@/components/viz/plotly'
import { MetaNode } from '@/spec/metanode'
import { ScoredDiseases, ScoredDrugs, ScoredGenes, ScoredPathways, ScoredPhenotypes, ScoredTissues } from '@/components/core/scored'
import { barchart_icon } from '@/icons'

export const BarChartsFromScoredT = [ScoredDiseases, ScoredDrugs, ScoredGenes, ScoredPathways, ScoredPhenotypes, ScoredTissues].map(ScoredT =>
  MetaNode(`BarChartFrom[${ScoredT.spec}]`)
    .meta({
      label: `Bar chart from ${ScoredT.meta.label}`,
      description: `Construct Bar Chart with ${ScoredT.meta.label}`,
      icon: [barchart_icon],
    })
    .inputs({ terms: ScoredT })
    .output(PlotlyPlot)
    .resolve(async (props) => await python(
      'components.data.barchart.createbarchart',
      { kargs: [props.inputs.terms], kwargs: { terms: ScoredT.meta.label } },
      message => props.notify({ type: 'info', message }),
    ))
    .story(props =>
      `To visualize the level of expression across ${ScoredT.meta.label.toLocaleLowerCase()}, a bar chart was created${''/* [FIGURE]*/}.`
    )
    .build()
)
