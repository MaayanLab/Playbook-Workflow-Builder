import { ScoredGenes } from '@/components/core/input/scored'
import { GeneSet } from '@/components/core/input/set'
import { idg_icon } from '@/icons'
import { MetaNode } from '@/spec/metanode'

export const IDGFilterT = [
  { kind: 'kinase', label: 'Kinase', },
  { kind: 'gpcr', label: 'GPCR', },
  { kind: 'ionchannel', label: 'Ion Channel', },
  { kind: 'darkkinase', label: 'Dark Kinase', },
  { kind: 'darkgpcr', label: 'Dark GPCR', },
  { kind: 'darkionchannel', label: 'Dark Ion Channel', },
].flatMap(({ kind, label }) => [
  MetaNode(`IDGFilter[${ScoredGenes.spec}, ${kind}]`)
    .meta({
      label: `Filter genes by ${label}`,
      description: `Filter genes by IDG ${label} classification`,
      icon: [idg_icon],
      pagerank: -1,
    })
    .inputs({ input: ScoredGenes })
    .output(ScoredGenes)
    .resolve(async (props) => {
      const req = await fetch('https://maayanlab.cloud/geneshot/api/idgclass')
      const { idg: { [kind]: filter } } = await req.json()
      const filterSet = new Set((filter as string[]).map(term => term.toLowerCase()))
      return props.inputs.input.filter(({ term }) => filterSet.has(term.toLocaleLowerCase()))
    })
    .story(props => `Genes were filtered by those classified as ${label} by IDG.`)
    .build(),
  MetaNode(`IDGFilter[${GeneSet.spec}, ${kind}]`)
    .meta({
      label: `Filter genes by ${label}`,
      description: `Filter genes by IDG ${label} classification`,
      icon: [idg_icon],
      pagerank: -1,
    })
    .inputs({ input: GeneSet })
    .output(GeneSet)
    .resolve(async (props) => {
      const req = await fetch('https://maayanlab.cloud/geneshot/api/idgclass')
      const { idg: { [kind]: filter } } = await req.json()
      const filterSet = new Set((filter as string[]).map(term => term.toLowerCase()))
      return {
        description: props.inputs.input.description ? `${props.inputs.input.description} filtered by ${label}` : undefined,
        set: props.inputs.input.set.filter((term) => filterSet.has(term.toLocaleLowerCase())),
      }
    })
    .story(props => `Genes were filtered by those classified as ${label} by IDG.`)
    .build(),
])
