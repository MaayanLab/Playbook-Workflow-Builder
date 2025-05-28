import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { ScoredGenes } from '@/components/core/scored'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import python from '@/utils/python'
import { archs4_icon, gtex_icon, targetranger_icon } from '@/icons'

export const TargtRangerScreenTargetsT = [
  { label: 'GTEx', bg: 'GTEx_transcriptomics', icon: gtex_icon, ref: '\\ref{doi:10.1038/ng.2653}', },
  { label: 'ARCHS4', bg: 'ARCHS4', icon: archs4_icon, ref: '\\ref{doi:10.1038/s41467-018-03751-6}', },
  ].map(({ label, bg, ref, icon }) =>
  MetaNode(`TargetRangerScreenTargets[${bg}]`)
    .meta({
      label: `Screen for Targets against ${label}`,
      description: `Identify significantly overexpressed genes when compared to normal tissue in ${label}`,
      icon: [targetranger_icon, icon],
      external: true,
    })
    .inputs({ input: GeneCountMatrix })
    .output(ScoredGenes)
    .resolve(async (props) => {
      return await python(
        'components.service.targetranger.targetscreener',
        { kargs: [props.inputs.input], kwargs: { bg } },
        message => props.notify({ type: 'info', message }),
      )
    })
    .story(props => ({
      abstract: `Significantly over-expressed genes when compared to tissue expression in ${label}${ref} were identified.`,
      introduction: `Several atlasing efforts profile human gene expression across tissues in both normal and diseased states. TargetRanger is a web server that compares uploaded RNA-seq expression data and identifies genes that are highly expressed when compared to various atlases\\ref{doi:10.1093/nar/gkad399}.`,
      methods: `The RNA-seq expression in ${props.input_refs?.input} is queried with TargetRanger\\ref{doi:10.1093/nar/gkad399} to identify highly expressed genes when compared with tissue expression in ${label}${ref}.`,
      tableLegend: `A table showing how significantly genes are highly expressed in ${props.input_refs?.input} when compared to tissue expression in ${label}${ref}.`,
    }))
    .build()
)
