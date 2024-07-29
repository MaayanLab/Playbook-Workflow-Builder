import python from '@/utils/python'
import { PlotlyPlot } from '@/components/viz/plotly'
import { MetaNode } from '@/spec/metanode'
import { ScoredDiseases, ScoredDrugs, ScoredGenes, ScoredPathways, ScoredPhenotypes, ScoredTissues } from '@/components/core/scored'
import { barplot_icon } from '@/icons'

export const BarChartsFromScoredT = [ScoredDiseases, ScoredDrugs, ScoredGenes, ScoredPathways, ScoredPhenotypes, ScoredTissues].map(ScoredT =>
  MetaNode(`BarChartFrom[${ScoredT.spec}]`)
    .meta({
      label: `Horizontal Bar Plot from ${ScoredT.meta.label}`,
      description: `Construct a horizontal bar plot with ${ScoredT.meta.label}`,
      icon: [barplot_icon],
    })
    .inputs({ terms: ScoredT })
    .output(PlotlyPlot)
    .resolve(async (props) => await python(
      'components.data.barchart.createbarchart',
      { kargs: [props.inputs.terms], kwargs: { terms: ScoredT.meta.label } },
      message => props.notify({ type: 'info', message }),
    ))
    .story(props => ({
      abstract: `To visualize the level of expression across ${ScoredT.meta.label.toLocaleLowerCase()}, a horizontal bar plot was created${''/* [FIGURE]*/}.`,
      legend: `A horizontal bar chart showing the gene expression levels across ${ScoredT.meta.label.toLowerCase()}.`,
    }))
    .build()
)
