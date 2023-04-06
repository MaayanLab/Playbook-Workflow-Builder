import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { ScoredGenes } from '@/components/core/input/scored'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import python from '@/utils/python'
import { targetranger_icon } from '@/icons'

export const TargtRangerScreenTargetsT = [
  { label: 'GTEx', bg: 'GTEx_transcriptomics' },
  { label: 'ARCHS4', bg: 'ARCHS4' },
].map(({ label, bg }) =>
  MetaNode(`TargetRangerScreenTargets[${bg}]`)
    .meta({
      label: `Screen for Targets against ${label}`,
      description: `Identify significantly overexpressed genes when compared to normal tissue in ${label}`,
      icon: [targetranger_icon],
    })
    .inputs({ input: GeneCountMatrix })
    .output(ScoredGenes)
    .resolve(async (props) => {
      return await python(
        'components.service.targetranger.targetscreener',
        { kargs: [props.inputs.input.url], kwargs: { bg } },
      )
    })
    .story(props => `Significantly over-expressed genes when compared to normal tissue in ${label} were identified.`)
    .build()
)
