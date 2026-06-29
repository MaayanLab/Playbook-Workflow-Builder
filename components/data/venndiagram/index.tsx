import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { norm_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'
import { GeneSet } from '@/components/core/set'


export const VennDiagram = MetaNode('VennDiagram')
  .meta({
    label: 'Venn Diagram from Gene Sets',
    description: 'Construct a ven diagram for 2-6 gene sets. The intersections of between each set will be displayed with a label for the size of each intersection.',
    icon: [norm_icon],
  })
  .inputs({sets:[GeneSet]})
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.venndiagram.createvenn',
    { kargs: [...props.inputs.sets]  },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `${props.inputs?.sets.length} sets were used to create a Venn Diagram.`,
    introduction: `Differential Gene Expression (DGE) anaysis is a technique used to identify gene that are expressed more or less between cells under two different experimental conditions. Analyzing the differences in gene expression levels can provide information about the role certain genes play in a cellular phenotype\\ref{doi:10.1186/gb-2013-14-9-r95}.`,
    methods: `Differential gene expression (DGE) analysis was performed between the Control group and Perturbation group samples. Based on the DGE analysis, each gene was assigned a log2-fold-change and statistical signifiance score. These scores were plotted using a volcano plot.`,
    legend: `A volcano plot showing the log2-fold-changes and statistical significance of each gene.`,
  }))
  .build()


  