import React from 'react'
import { MetaNode } from '@/spec/metanode'

const Plot = React.lazy(() => import('./bokeh'))

export type BokehJson = any

export const BokehPlot = MetaNode.createData('BokehPlot')
  .meta({
    label: 'Bokeh Plot',
    description: 'A plot rendered using the bokeh library',
  })
  .codec<BokehJson>()
  .view(props => <Plot plot={props} />)
  .build()
