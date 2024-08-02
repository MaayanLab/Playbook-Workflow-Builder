import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { norm_icon } from '@/icons'


export const QuantileNormalizeGeneCountMatrix = MetaNode('QuantileNormalizeGeneCountMatrix')
  .meta({
    label: 'Quantile Normalize A Gene Count Matrix',
    description: 'Quantile normalize a gene count matrix, return a gene count matrix',
    icon: [norm_icon],
  })
  .inputs({ matrix: GeneCountMatrix })
  .output(GeneCountMatrix)
  .resolve(async (props) => await python(
    'components.data.quantile_normalization.quantile_normalize_gene_count_matrix',
    { kargs: [props.inputs.matrix]  },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The gene count matrix was then quantile normalized\\ref{doi:10.1038/s41598-020-72664-6}.`,
    introduction: `Quantile normalization (QN) is technique often used with high-dimensional -omics data including RNA-sequencing and proteomics. It assumes that sample data follow similar distributions regardless of the sample class. Class-specific QN has been demonstrated to optimize performance compared to whole-data QN\\ref{doi:10.1038/s41598-020-72664-6}.`,
    methods: `The gene count matrix is then normalized using class-specific quantile normalization\\ref{doi:10.1038/s41598-020-72664-6}.`,
    legend: `A table displaying the quantile-normalized gene count matrix.`,
  }))
  .build()
