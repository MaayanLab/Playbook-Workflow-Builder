import { MetaNode } from '@/spec/metanode'
import dynamic from 'next/dynamic'
import { z } from 'zod'
import { plot_icon } from '@/icons'

const CytoscapeCanvas = dynamic(() => import('./cytoscape'), { ssr: false, loading: () => <div className="prose max-w-none">Loading...</div> })

export const GraphPlot = MetaNode('GraphPlot')
  .meta({
    label: 'Graph Plot',
    description: 'A graph plot rendered using the [cytoscape.js library](https://js.cytoscape.org/)',
    icon: [plot_icon],
  })
  .codec(z.object({
    nodes: z.array(z.object({
      id: z.string(),
      label: z.string().optional(),
      type: z.string(),
      color: z.string().optional(),
    })),
    edges: z.array(z.object({
      source: z.string(),
      target: z.string(),
    })),
  }))
  .view(props => <CytoscapeCanvas {...props} />)
  .build()
