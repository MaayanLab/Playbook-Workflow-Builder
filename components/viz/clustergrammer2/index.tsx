import React from 'react'
import { MetaNode } from '@/spec/metanode'
import dynamic from 'next/dynamic'
import { plot_icon } from '@/icons'
import type { Clustergrammer2Data } from './clustergrammer2'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import python from '@/utils/python'
import { GMT } from '@/components/data/gene_matrix_transpose'

const Plot = dynamic(() => import('./clustergrammer2'), { ssr: false, loading: () => <div>Loading...</div> })

export const Clustergrammer2Plot = MetaNode('Clustergrammer2Plot')
  .meta({
    label: 'Clustergrammer2 Plot',
    description: 'A plot rendered using the clustergrammer-gl library',
    icon: [plot_icon],
  })
  .codec<Clustergrammer2Data>()
  .view(props => <Plot {...props} />)
  .build()

export const Clustergrammer2FromGeneCountMatrix = MetaNode('Clustergrammer2FromGeneCountMatrix')
  .meta({
    label: 'Construct Clustergrammer2 Plot',
    description: 'Plot matrix using clustergrammer2',
  })
  .inputs({ data: GeneCountMatrix })
  .output(Clustergrammer2Plot)
  .resolve(async (props) => await python(
    'components.viz.clustergrammer2.clustergrammer2_from_df',
    { kargs: [props.inputs.data.url] },
  ))
  .story(props => `Clustergrammer2 was used to visualize the data.`)
  .build()


export const Clustergrammer2FromGeneMatrixTranspose = MetaNode('Clustergrammer2FromGeneMatrixTranspose')
  .meta({
    label: 'Construct Clustergrammer2 Plot',
    description: 'Plot matrix using clustergrammer2',
  })
  .inputs({ data: GMT })
  .output(Clustergrammer2Plot)
  .resolve(async (props) => await python(
    'components.viz.clustergrammer2.clustergrammer2_from_gmt',
    { kargs: [props.inputs.data] },
  ))
  .story(props => `Clustergrammer2 was used to visualize the data.`)
  .build()
