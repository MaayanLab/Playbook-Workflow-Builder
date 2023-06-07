import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { DrugSet, GeneSet } from '@/components/core/input/set'
import { Disease, Drug, Gene, Pathway, Phenotype, Tissue } from '@/components/core/primitives'
import { DiseaseTerm, DrugTerm, GeneTerm, PathwayTerm, PhenotypeTerm, TissueTerm } from '@/components/core/input/term'
import { ScoredDrugs, ScoredGenes, ScoredDiseases, ScoredPathways, ScoredPhenotypes, ScoredTissues } from '@/components/core/input/scored'
import { Table, Cell, Column } from '@/app/components/Table'
import classNames from 'classnames'

function withPrecision(value: number | string, precision: number) {
  if (typeof value === 'number') return value.toPrecision(precision)
  else return value
}

export const TopKScoredT = [
  { ScoredT: ScoredDrugs, TermT: DrugTerm, T: Drug },
  { ScoredT: ScoredGenes, TermT: GeneTerm, T: Gene },
  { ScoredT: ScoredDiseases, TermT: DiseaseTerm, T: Disease },
  { ScoredT: ScoredPathways, TermT: PathwayTerm, T: Pathway },
  { ScoredT: ScoredPhenotypes, TermT: PhenotypeTerm, T: Phenotype },
  { ScoredT: ScoredTissues, TermT: TissueTerm, T: Tissue },
].flatMap(({ ScoredT, TermT, T }) => [
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
    .build(),
  MetaNode(`OneScoredT[${ScoredT.spec}]`)
    .meta({
      label: `Select One ${TermT.meta.label}`,
      description: `Select one ${TermT.meta.label}`,
    })
    .inputs({ scored: ScoredT })
    .output(TermT)
    .prompt(props => {
      const scored = props.inputs.scored
      const [selected, setSelected] = React.useState<number | undefined>()
      React.useEffect(() => {
        if (props.output !== undefined) {
          const index = scored.findIndex(({ term }) => term === props.output)
          if (index >= 0) {
            setSelected(index)
            return
          }
        }
        if (scored.length > 0) {
          setSelected(0)
        }
      }, [props.output, scored])
      return (
        <div>
          <Table
            height={500}
            cellRendererDependencies={[scored]}
            rowHeaderCellRenderer={(row) =>
              <div
                className="text-center block"
                onClick={evt => {
                  setSelected(row)
                  props.submit(scored[row].term)
                }}
              >
                <input type="radio" checked={selected === row} />
              </div>
            }
            numRows={scored.length}
            enableGhostCells
            enableFocusedCell
          >
            <Column
              name={T.label}
              cellRenderer={row => <Cell key={row+''}>{scored[row].term}</Cell>}
            />
            <Column
              name="ZScore"
              cellRenderer={row => <Cell key={row+''}>{withPrecision(scored[row].zscore, 3)}</Cell>}
            />
          </Table>
        </div>
      )
    })
    .story(props => props.output ? `${props.output} was chosen for further investigation.` : '')
    .build()
])

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
