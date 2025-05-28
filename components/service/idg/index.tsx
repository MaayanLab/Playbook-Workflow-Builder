import { ScoredGenes } from '@/components/core/scored'
import { GeneSet } from '@/components/core/set'
import { idg_icon } from '@/icons'
import { MetaNode } from '@/spec/metanode'
import * as dict from '@/utils/dict'

async function idg_filter(kind: string) {
  const req = await fetch('https://maayanlab.cloud/geneshot/api/idgclass')
  let filter: string[]
  const { idg } = await req.json()
  if (kind === 'all') {
    filter = dict.values(
      dict.filter(idg, ({ key }) => (key as string).startsWith('dark'))
    ).reduce((all, proteins) => [...all, ...proteins], [])
  } else {
    filter = idg[kind]
  }
  return new Set((filter || [] as string[]).map(term => term.toLowerCase()))
}

export const IDGFilterT = [
  { kind: 'all', label: 'Understudied Proteins', },
  { kind: 'kinase', label: 'Kinases', },
  { kind: 'gpcr', label: 'GPCRs', },
  { kind: 'ionchannel', label: 'Ion Channels', },
  { kind: 'darkkinase', label: 'Dark Kinases', },
  { kind: 'darkgpcr', label: 'Dark GPCRs', },
  { kind: 'darkionchannel', label: 'Dark Ion Channels', },
].flatMap(({ kind, label }) => [
  MetaNode(`IDGFilter[${ScoredGenes.spec}, ${kind}]`)
    .meta({
      label: `Filter genes by ${label}`,
      description: `Based on IDG proteins list`,
      icon: [idg_icon],
      pagerank: -1,
      external: true,
    })
    .inputs({ input: ScoredGenes })
    .output(ScoredGenes)
    .resolve(async (props) => {
      const filterSet = await idg_filter(kind)
      return props.inputs.input.filter(({ term }) => filterSet.has(term.toLocaleLowerCase()))
    })
    .story(props => ({
      abstract: `Genes were filtered by IDG ${label}\\ref{IDG Protein List, https://druggablegenome.net/IDGProteinList}.`,
      introduction: `The Illuminating the Druggable Genome (IDG) Program seeks to improve our understanding of the properties and functions of proteins not currently well studied within commonly drug-targeted protein families. IDG maintains a list of understudied proteins from three key druggable protein families (GPCR, Ion Channels and Kinases).`,
      methods: `Genes in ${props.input_refs?.input} are filtered by the list of IDG ${label}.`,
      tableLegend: `A table of IDG ${label} filtered from ${props.input_refs?.input}.`,
    }))
    .build(),
  MetaNode(`IDGFilter[${GeneSet.spec}, ${kind}]`)
    .meta({
      label: `Filter genes by IDG ${label}`,
      description: `Based on IDG understudied proteins list`,
      icon: [idg_icon],
      pagerank: -1,
      external: true,
    })
    .inputs({ input: GeneSet })
    .output(GeneSet)
    .resolve(async (props) => {
      const filterSet = await idg_filter(kind)
      return {
        description: props.inputs.input.description ? `${props.inputs.input.description} filtered by ${label}` : undefined,
        set: props.inputs.input.set.filter((term) => filterSet.has(term.toLocaleLowerCase())),
      }
    })
    .story(props => ({
      abstract: `Genes were filtered by IDG ${label}\\ref{IDG Understudied Proteins, https://druggablegenome.net/AboutIDGProteinList}.`,
      introduction: `The Illuminating the Druggable Genome (IDG) Program seeks to improve our understanding of the properties and functions of proteins not currently well studied within commonly drug-targeted protein families. IDG maintains a list of understudied proteins from three key druggable protein families (GPCR, Ion Channels and Kinases).`,
      methods: `Genes in ${props.input_refs?.input} are filtered by the list of IDG ${label}.`,
      tableLegend: `A table of IDG ${label} filtered from ${props.input_refs?.input}.`,
    }))
    .build(),
])
