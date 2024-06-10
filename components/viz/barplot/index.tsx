import python from '@/utils/python'
import { PlotlyPlot } from '@/components/viz/plotly'
import { MetaNode } from '@/spec/metanode'
import { ScoredDiseases, ScoredDrugs, ScoredGenes, ScoredPathways, ScoredPhenotypes, ScoredTissues } from '@/components/core/scored'
import { barplot_icon } from '@/icons'

export const BarplotFromScoredT = [ScoredDiseases, ScoredDrugs, ScoredGenes, ScoredPathways, ScoredPhenotypes, ScoredTissues].map(ScoredT =>
  MetaNode(`BarplotFrom[${ScoredT.spec}]`)
    .meta({
      label: `Vertical bar plot from ${ScoredT.meta.label}`,
      description: `Construct a vertical bar plot with ${ScoredT.meta.label}`,
      icon: [barplot_icon],
    })
    .inputs({ terms: ScoredT })
    .output(PlotlyPlot)
    .resolve(async (props) => await python(
      'components.viz.barplot.barplot',
      { kargs: [props.inputs.terms], kwargs: { terms: ScoredT.meta.label } },
      message => props.notify({ type: 'info', message }),
    ))
    .story(props => ({
      abstract: `To visualize the level of expression across ${ScoredT.meta.label.toLocaleLowerCase()}, a vertical bar plot was created${''/* [FIGURE]*/}.`
    }))
    .build()
)
