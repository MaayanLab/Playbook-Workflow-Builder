import { ScoredGenes } from '@/components/core/scored'
import { GeneSet } from '@/components/core/set'
import { idg_icon } from '@/icons'
import { MetaNode } from '@/spec/metanode'
import Citable from '@/utils/citations'
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
      description: `Based on IDG understudied proteins list`,
      icon: [idg_icon],
      pagerank: -1,
    })
    .inputs({ input: ScoredGenes })
    .output(ScoredGenes)
    .resolve(async (props) => {
      const filterSet = await idg_filter(kind)
      return props.inputs.input.filter(({ term }) => filterSet.has(term.toLocaleLowerCase()))
    })
    .story(props => ({ abstract: Citable.text`Genes were filtered by IDG ${label} [${Citable.cite('IDG Understudied Proteins, https://druggablegenome.net/AboutIDGProteinList')}].` }))
    .build(),
  MetaNode(`IDGFilter[${GeneSet.spec}, ${kind}]`)
    .meta({
      label: `Filter genes by IDG ${label}`,
      description: `Based on IDG understudied proteins list`,
      icon: [idg_icon],
      pagerank: -1,
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
    .story(props => ({ abstract: Citable.text`Genes were filtered by IDG ${label} [${Citable.cite('IDG Understudied Proteins, https://druggablegenome.net/AboutIDGProteinList')}].` }))
    .build(),
])
