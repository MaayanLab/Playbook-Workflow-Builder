import React from 'react'
import { MetaNode } from '@/spec/metanode'
import type { PlotParams } from 'react-plotly.js'

const Plot = React.lazy(() => import('react-plotly.js'))

export type PlotlyJson = {
  data: PlotParams['data'],
  layout: PlotParams['layout'],
  frames?: PlotParams['frames'],
}

export const PlotlyPlot = MetaNode.createData('PlotlyPlot')
  .meta({
    label: 'Plotly Plot',
    description: 'A plot rendered using the plotly.js library',
  })
  .codec<PlotlyJson>()
  .view(props => <Plot {...props} />)
  .build()
