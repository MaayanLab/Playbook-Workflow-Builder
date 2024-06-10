import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { MetaboliteCountMatrix } from '@/components/data/metabolite_count_matrix'
import { norm_icon } from '@/icons'


export const ZScoreNormalizeMetaboliteCountMatrix = MetaNode('ZScoreNormalizeMetaboliteCountMatrix')
  .meta({
    label: 'Z Score Normalize A Metabolite Count Matrix',
    description: 'Z-score normalize a metabolite count matrix, return a metabolite count matrix',
    icon: [norm_icon],
  })
  .inputs({ matrix: MetaboliteCountMatrix })
  .output(MetaboliteCountMatrix)
  .resolve(async (props) => await python(
    'components.data.met_z_score_normalization.z_score_normalize_metabolite_count_matrix',
    { kargs: [props.inputs.matrix]  },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The metabolite count matrix was then Z-score normalized.`,
  }))
  .build()
