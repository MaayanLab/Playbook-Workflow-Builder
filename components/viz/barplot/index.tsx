import python from '@/utils/python'
import { PlotlyPlot } from '@/components/viz/plotly'
import { MetaNode } from '@/spec/metanode'
import { ScoredDiseases, ScoredDrugs, ScoredGenes, ScoredPathways, ScoredPhenotypes, ScoredTissues } from '@/components/core/input/scored'
import { barchart_icon } from '@/icons'

export const BarplotFromScoredT = [ScoredDiseases, ScoredDrugs, ScoredGenes, ScoredPathways, ScoredPhenotypes, ScoredTissues].map(ScoredT =>
  MetaNode(`BarplotFrom[${ScoredT.spec}]`)
    .meta({
      label: `Barplot from ${ScoredT.meta.label}`,
      description: `Construct Barplot with ${ScoredT.meta.label}`,
      icon: [barchart_icon],
    })
    .inputs({ terms: ScoredT })
    .output(PlotlyPlot)
    .resolve(async (props) => await python(
      'components.viz.barplot.barplot',
      { kargs: [props.inputs.terms], kwargs: { terms: ScoredT.meta.label } },
    ))
    .story(props =>
      `To visualize the level of expression across ${ScoredT.meta.label.toLocaleLowerCase()}, a bar plot was created${''/* [FIGURE]*/}.`
    )
    .build()
)
