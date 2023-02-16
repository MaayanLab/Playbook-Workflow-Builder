import React from 'react'
import { MetaNode } from '@/spec/metanode'
import dynamic from 'next/dynamic'

const Plot = dynamic(() => import('./bokeh'), { ssr: false })

export type BokehJson = any

export const BokehPlot = MetaNode('BokehPlot')
  .meta({
    label: 'Bokeh Plot',
    description: 'A plot rendered using the bokeh library',
  })
  .codec<BokehJson>()
  .view(props => <Plot plot={props} />)
  .build()
