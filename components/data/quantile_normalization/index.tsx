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
  ))
  .story(props =>
    `The gene count matrix was then quantile normalized [\\ref{doi:10.1038/s41598-020-72664-6}].`
  )
  .build()
