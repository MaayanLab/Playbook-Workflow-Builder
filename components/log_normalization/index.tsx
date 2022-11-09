import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneCountMatrix } from '@/components/gene_count_matrix'


export const LogNormalizeGeneCountMatrix = MetaNode.createProcess('LogNormalizeGeneCountMatrix')
  .meta({
    label: 'Log Normalize A Gene Count Matrix',
    description: 'Log_2 normalize a gene count matrix, return a gene count matrix',
  })
  .codec()
  .inputs({ matrix: GeneCountMatrix })
  .output(GeneCountMatrix)
  .resolve(async (props) => await python(
    'components.log_normalization.log_normalize_gene_count_matrix',
    { kargs: [props.inputs.matrix]  },
  ))
  .build()
