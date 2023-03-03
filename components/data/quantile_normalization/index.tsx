import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'


export const QuantileNormalizeGeneCountMatrix = MetaNode('QuantileNormalizeGeneCountMatrix')
  .meta({
    label: 'Quantile Normalize A Gene Count Matrix',
    description: 'Quantile normalize a gene count matrix, return a gene count matrix',
  })
  .inputs({ matrix: GeneCountMatrix })
  .output(GeneCountMatrix)
  .resolve(async (props) => await python(
    'components.data.quantile_normalization.quantile_normalize_gene_count_matrix',
    { kargs: [props.inputs.matrix]  },
  ))
  .story(props =>
    `The gene count matrix was then quantile normalized.`
  )
  .build()
