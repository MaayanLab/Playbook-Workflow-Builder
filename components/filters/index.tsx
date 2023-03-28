import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { DrugSet, GeneSet } from '@/components/core/input/set'
import { ScoredDrugs, ScoredGenes, ScoredDiseases, ScoredPathways, ScoredPhenotypes, ScoredTissues } from '@/components/core/input/scored'
import classNames from 'classnames'

export const TopKScoredT = [
  ScoredDrugs,
  ScoredGenes,
  ScoredDiseases,
  ScoredPathways,
  ScoredPhenotypes,
  ScoredTissues,
].map((ScoredT) =>
  MetaNode(`TopKScoredT[${ScoredT.spec}]`)
    .meta({
      label: `Top ${ScoredT.meta.label}`,
      description: `Select the top ${ScoredT.meta.label}`,
    })
    .inputs({ scored: ScoredT })
    .output(ScoredT)
    .prompt(props => {
      const topK = React.useCallback((scored: typeof props.inputs.scored, k: number) => {
        scored.sort((a, b) =>
          a.zscore === b.zscore ? 0
          : a.zscore === 'inf' ? 1
          : b.zscore === 'inf' ? -1
          : a.zscore === '-inf' || a.zscore === 'nan' ? -1
          : b.zscore === '-inf' || b.zscore === 'nan' ? 1
          : a.zscore < b.zscore ? 1
          : -1
        )
        return scored.slice(0, k)
      }, [])
      return (
        <div className="flex flex-col gap-2">
          <div className="flex flex-row gap-2">
            {[10, 20, 50].map(k =>
              <button
                key={k}
                className={classNames('btn', {'btn-primary': props.output?.length === k, 'btn-secondary': props.output?.length !== k})}
                onClick={() => {props.submit(topK(props.inputs.scored, k))}}
              >Top {k}</button>
            )}
          </div>
          {ScoredT.view(props.inputs.scored)}
        </div>
      )
    })
    .story(props => `The top ${props.output?.length || 'K'} ${ScoredT.meta.label} were selected.`)
    .build()
)

export const SetFromScoredT = [
  { ScoredT: ScoredDrugs, SetT: DrugSet, },
  { ScoredT: ScoredGenes, SetT: GeneSet, },
].map(({ ScoredT, SetT }) =>
  MetaNode(`SetFromScored[${ScoredT.spec}]`)
    .meta({
      label: `Set from ${ScoredT.meta.label}`,
      description: `Construct Set with ${ScoredT.meta.label}`,
    })
    .inputs({ scored: ScoredT })
    .output(SetT)
    .resolve(async (props) => ({ set: props.inputs.scored.map(({ term }) => term) }))
    .story(props =>
      `A set was constructed using the ${ScoredT.meta.label.toLocaleLowerCase()}.`
    )
    .build()
)
