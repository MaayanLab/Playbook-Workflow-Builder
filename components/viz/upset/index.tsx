import { MetaNode } from '@/spec/metanode'
import { plot_icon } from '@/icons'
import { GMT } from '@/components/data/gene_matrix_transpose'
import { DMT } from '@/components/data/drug_matrix_transpose'
import * as dict from '@/utils/dict'
import dynamic from 'next/dynamic'

const UpsetPlotV2 = dynamic(() => import('./UpSet').then(({ UpsetPlotV2 }) => UpsetPlotV2), { ssr: false })

export const UpSetPlot = MetaNode(`UpSetPlot`)
  .meta({
    label: `UpSet plot`,
    description: 'An UpSet plot',
    icon: [plot_icon],
  })
  .codec(GMT.codec)
  .view(xmt => <UpsetPlotV2
      selectedSets={dict.items(xmt).map(({ key, value }) => ({ alphabet: key, genes: value.set.map(gene_symbol => ({ gene_symbol }))}))}
      setOverlap={() => {}}
    />)
  .build()

export const UpSetFromXMT = [GMT, DMT].map(XMT =>
  MetaNode(`UpSetFrom[${XMT.spec}]`)
    .meta({
      label: `UpSet plot from ${XMT.meta.label}`,
      description: `Construct UpSet plot with ${XMT.meta.label}`,
      icon: [plot_icon],
    })
    .inputs({ matrix: XMT })
    .output(UpSetPlot)
    .resolve(async (props) => props.inputs.matrix)
    .story(props => ({
      abstract: `A UpSet plot was constructed with the ${XMT.meta.label.toLocaleLowerCase()}.`,
      introduction: `UpSet is a visualization technique for the quantitative analysis of sets, their intersections, and aggregates of intersections\\ref{doi:10.1109/TVCG.2014.2346248}.`,
      methods: `The sets in the gene matrix transpose are visualized using a D3-based interactive UpSet plot\\ref{doi:10.1109/TVCG.2014.2346248}\\ref{d3 UpSet Plot, https://github.com/chuntul/d3-upset}.`,
      legend: `An interactive UpSet\\ref{doi:10.1109/TVCG.2014.2346248} plot visualizing the ${XMT.meta.label.toLocaleLowerCase()} in${''/* {props.input_ref.matrix}*/}.`,
    }))
    .build()
)
