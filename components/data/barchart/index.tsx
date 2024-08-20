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
      abstract: `To visualize the ${ScoredT.meta.label.toLocaleLowerCase()}, a horizontal bar plot was created (${props.output_ref}).`,
      introduction: `Plotly.js is a browser-based graphing library for Python which facilitates the construction of interactive vector graphic charts \\ref{Plotly Technologies Inc. Collaborative data science. Montreal, QC, 2015. https://plot.ly}  `,
      methods: `The table in ${props.input_refs?.terms} is visualized with a horizontal bar graph using Plotly.js\\ref{Plotly Technologies Inc. Collaborative data science. Montreal, QC, 2015. https://plot.ly}.`,
      legend: `A horizontal bar plot visualizing the values from ${props.input_refs?.terms}.`,
    }))
    .build()
)
