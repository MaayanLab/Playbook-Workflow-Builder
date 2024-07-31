import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneSignature } from '@/components/data/gene_signature'
import { norm_icon } from '@/icons'
import { PlotlyPlot } from '@/components/viz/plotly'



export const VolcanoPlot = MetaNode('VolcanoPlot')
  .meta({
    label: 'Volcano Plot from Differential Expression Table',
    description: 'Construct a scatter plot which displays the log2-fold changes and statistical significance of each gene calculated by performing differential gene expression analysis comparing samples in the Control group to samples in the Perturbation group. Every point in the plot represents a gene.',
    icon: [norm_icon],
  })
  .inputs({ sig: GeneSignature })
  .output(PlotlyPlot)
  .resolve(async (props) => await python(
    'components.data.volcanoplot.createvolcano',
    { kargs: [props.inputs.sig]  },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The data in the differential expression table was then visualized as a volcano plot.`,
    introduction: `Differential Gene Expression (DGE) anaysis is a technique used to identify gene that are expressed more or less between cells under two different experimental conditions. Analyzing the differences in gene expression levels can provide information about the role certain genes play in a cellular phenotype\\ref{doi:10.1186/gb-2013-14-9-r95}.`,
    methods: `Differential gene expression (DGE) analysis was performed between the Control group and Perturbation group samples. Based on the DGE analysis, each gene was assigned a log2-fold-change and statistical signifiance score. These scores were plotted using a volcano plot.`,
    legend: `A volcano plot showing the log2-fold-changes and statistical significance of each gene.`,
  }))
  .build()


  