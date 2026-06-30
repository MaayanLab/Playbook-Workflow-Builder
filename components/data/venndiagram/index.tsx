import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { supervenn_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'
import { GeneSet } from '@/components/core/set'
import { GMT } from '../gene_matrix_transpose'

export const VennDiagram = MetaNode('VennDiagram')
  .meta({
    label: 'Venn Diagram from Gene Sets',
    description: 'Construct a ven diagram for 2-6 gene sets. The intersections of between each set will be displayed with a label for the size of each intersection.',
    icon: [supervenn_icon],
  })
  .inputs({sets:[GeneSet]})
  .output(PlotlyPlot)
  .resolve(async (props) => {
    if ((props.inputs.sets.length < 2) || (props.inputs.sets.length >6)) {
      throw new Error("Please submit between 2 and 6 gene sets.")
    }
    
    return await python(
    'components.data.venndiagram.createvenn',
    { kargs: [...props.inputs.sets]  },
    message => props.notify({ type: 'info', message }),
  )})
  .story(props => ({
    abstract: `${props.inputs?.sets.length} sets were used to create a Venn Diagram.`,
    legend: `A venn diagram of input gene sets.`,
  }))
  .build()

export const VennDiagramGMT = MetaNode('VennDiagram')
  .meta({
    label: 'Venn Diagram from Gene Set Library',
    description: 'Construct a ven diagram for 2-6 gene sets. The intersections of between each set will be displayed with a label for the size of each intersection.',
    icon: [supervenn_icon],
  })
  .inputs({gmt:GMT})
  .output(PlotlyPlot)
  .resolve(async (props) => {
    const sets = Object.entries(props.inputs.gmt).map(([term, { set }]) => ({set,description: term,}))

    if ((sets.length < 2) || (sets.length >6)) {
      throw new Error("Please submit between 2 and 6 gene sets.")
    }

    return await python(
      'components.data.venndiagram.createvenn',
      { kargs: [...sets]  },
      message => props.notify({ type: 'info', message }),
    )
  })
  .story(props => ({
    abstract: `${props.inputs ? Object.keys(props.inputs.gmt).length + " sets": "Sets"} were used to create a Venn Diagram.`,
    legend: `A venn diagram of input gene sets.`,
  }))
  .build()
  