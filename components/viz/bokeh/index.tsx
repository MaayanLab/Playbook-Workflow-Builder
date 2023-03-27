import React from 'react'
import { MetaNode } from '@/spec/metanode'
import dynamic from 'next/dynamic'
import { plot_icon } from '@/icons'

const Plot = dynamic(() => import('./bokeh'), { ssr: false })

export type BokehJson = any

export const BokehPlot = MetaNode('BokehPlot')
  .meta({
    label: 'Bokeh Plot',
    description: 'A figure created with the [Bokeh Library](https://bokeh.org/)',
    icon: [plot_icon],
  })
  .codec<BokehJson>()
  .view(props => <Plot plot={props} />)
  .build()
