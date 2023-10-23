import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { norm_icon } from '@/icons'


export const LogNormalizeGeneCountMatrix = MetaNode('LogNormalizeGeneCountMatrix')
  .meta({
    label: 'Log Normalize A Gene Count Matrix',
    description: 'Log_2 normalize a gene count matrix, return a gene count matrix',
    icon: [norm_icon],
  })
  .inputs({ matrix: GeneCountMatrix })
  .output(GeneCountMatrix)
  .resolve(async (props) => await python(
    'components.data.log_normalization.log_normalize_gene_count_matrix',
    { kargs: [props.inputs.matrix]  },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props =>
    `The gene count matrix was then log2 transformed.`
  )
  .build()
