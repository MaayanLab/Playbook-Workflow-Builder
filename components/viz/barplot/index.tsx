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
      abstract: `To visualize the ${ScoredT.meta.label.toLocaleLowerCase()}, a vertical bar plot was created ${props.output_ref}.`,
      introduction: `Plotly.js is a browser-based graphing library for Python which facilitates the construction of interactive vector graphic charts \\ref{Plotly Technologies Inc. Collaborative data science. Montréal, QC, 2015. https://plot.ly}`,
      methods: `The table in ${props.input_refs?.terms} are visualized with a vertical bar graph using Plotly.js\\ref{Plotly Technologies Inc. Collaborative data science. Montréal, QC, 2015. https://plot.ly}.`,
      legend: `A vertical bar plot visualizing the values from ${props.input_refs?.terms}.`,
    }))
    .build()
)
