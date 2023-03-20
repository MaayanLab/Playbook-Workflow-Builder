import { MetaNode } from '@/spec/metanode'
import { DrugSet, GeneSet } from '@/components/core/input/set'
import { ScoredDrugs, ScoredGenes } from '@/components/core/input/scored'
  
export const SetFromScoredT = [
  { ScoredT: ScoredDrugs, SetT: DrugSet, },
  { ScoredT: ScoredGenes, SetT: GeneSet, },
].map(({ ScoredT, SetT }) =>
  MetaNode(`SetFromTop50[${ScoredT.spec}]`)
    .meta({
      label: `Set from ${ScoredT.meta.label}`,
      description: `Construct Set with top 50 ${ScoredT.meta.label}`,
    })
    .inputs({ scored: ScoredT })
    .output(SetT)
    .resolve(async (props) => {
      const scored = props.inputs.scored
      scored.sort((a, b) =>
        a.zscore === b.zscore ? 0
        : a.zscore === 'inf' ? 1
        : b.zscore === 'inf' ? -1
        : a.zscore === '-inf' || a.zscore === 'nan' ? -1
        : b.zscore === '-inf' || b.zscore === 'nan' ? 1
        : a.zscore < b.zscore ? 1
        : -1
      )
      return { set: scored.slice(0, 50).map(({ term }) => term) }
    })
    .story(props =>
      `A set was constructed using the top 50 ${ScoredT.meta.label.toLocaleLowerCase()}.`
    )
    .build()
)
