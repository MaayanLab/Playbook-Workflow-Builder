import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { norm_icon } from '@/icons'


export const ZScoreNormalizeGeneCountMatrix = MetaNode('ZScoreNormalizeGeneCountMatrix')
  .meta({
    label: 'Z Score Normalize A Gene Count Matrix',
    description: 'Z-score normalize a gene count matrix, return a gene count matrix',
    icon: [norm_icon],
  })
  .inputs({ matrix: GeneCountMatrix })
  .output(GeneCountMatrix)
  .resolve(async (props) => await python(
    'components.data.z_score_normalization.z_score_normalize_gene_count_matrix',
    { kargs: [props.inputs.matrix]  },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The gene count matrix was then Z-score normalized.`,
    introduction: `Gene count are normalized to control for different levels of variance depending on the gene expression levels.`,
    methods: `The gene count matrix was normalized using a z-score transformation.`,
    tableLegend: `A table displaying the z-score normalized gene count matrix.`,
  }))
  .build()
