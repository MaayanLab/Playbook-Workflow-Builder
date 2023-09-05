import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { DrugSet, GeneSet } from '@/components/core/input/set'
import { Disease, Drug, Gene, Pathway, Phenotype, Tissue } from '@/components/core/primitives'
import { DiseaseTerm, DrugTerm, GeneTerm, PathwayTerm, PhenotypeTerm, TissueTerm } from '@/components/core/input/term'
import { ScoredDrugs, ScoredGenes, ScoredDiseases, ScoredPathways, ScoredPhenotypes, ScoredTissues } from '@/components/core/input/scored'
import { Table, Cell, Column } from '@/app/components/Table'
import * as dict from '@/utils/dict'
import * as math from '@/utils/math'
import classNames from 'classnames'
import pluralize from 'pluralize'

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
  { ScoredT: ScoredDrugs, SetT: DrugSet, TermT: DrugTerm, T: Drug, },
  { ScoredT: ScoredGenes, SetT: GeneSet, TermT: GeneTerm, T: Gene, },
].flatMap(({ ScoredT, SetT, TermT, T }) => [
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
    .build(),
  MetaNode(`OneSetT[${SetT.spec}]`)
    .meta({
      label: `Select One ${TermT.meta.label}`,
      description: `Select one ${TermT.meta.label}`,
    })
    .inputs({ set: SetT })
    .output(TermT)
    .prompt(props => {
      const set = props.inputs.set.set
      const [selected, setSelected] = React.useState<number | undefined>()
      React.useEffect(() => {
        if (props.output !== undefined) {
          const index = set.findIndex((term) => term === props.output)
          if (index >= 0) {
            setSelected(index)
            return
          }
        }
        if (set.length > 0) {
          setSelected(0)
        }
      }, [props.output, set])
      return (
        <div>
          <Table
            height={500}
            cellRendererDependencies={[set]}
            rowHeaderCellRenderer={(row) =>
              <div
                className="text-center block"
                onClick={evt => {
                  setSelected(row)
                  props.submit(set[row])
                }}
              >
                <input type="radio" checked={selected === row} />
              </div>
            }
            numRows={set.length}
            enableGhostCells
            enableFocusedCell
          >
            <Column
              name={T.label}
              cellRenderer={row => <Cell key={row+''}>{set[row]}</Cell>}
            />
          </Table>
        </div>
      )
    })
    .story(props => props.output ? `${props.output} was chosen for further investigation.` : '')
    .build(),
  MetaNode(`SomeSetT[${SetT.spec}]`)
    .meta({
      label: `Select Some ${pluralize(TermT.meta.label)}`,
      description: `Select some ${pluralize(TermT.meta.label)}`,
    })
    .inputs({ set: SetT })
    .output(SetT)
    .prompt(props => {
      const set = props.inputs.set.set
      const [selected, setSelected] = React.useState({} as Record<string, true>)
      React.useEffect(() => {
        if (props.output !== undefined) {
          const selected_ = {} as Record<string, true>
          set.forEach(item => {
            if (!props.output?.set.includes(item)) {
              selected_[item] = true
            }
          })
          setSelected(selected_)
        }
      }, [props.output, set])
      return (
        <div>
          <Table
            height={500}
            cellRendererDependencies={[set]}
            rowHeaderCellRenderer={(row) =>
              <div
                className="text-center block"
                onClick={evt => {
                  setSelected(({ [set[row]]: current, ...selected }) => !current ? ({ ...selected, [set[row]]: true }) : selected)
                }}
              >
                <input type="checkbox" checked={!selected[set[row]]} />
              </div>
            }
            numRows={set.length}
            shape={[set.length - Object.keys(selected).length]}
            enableGhostCells
            enableFocusedCell
          >
            <Column
              name={T.label}
              cellRenderer={row => <Cell key={row+''}>{set[row]}</Cell>}
            />
          </Table>
          <button className="bp4-button bp4-large" onClick={async () => {
            props.submit({
              description: props.inputs.set.description ? `Filtered ${props.inputs.set.description}` : undefined,
              set: props.inputs.set.set.filter((item) => !selected[item])
            })
          }}>Submit</button>
        </div>
      )
    })
    .story(props => props.output ? `Some genes were selected for further investigation.` : '')
    .build(),
])

export const ReduceMultiScoredT = [
  { ScoredT: ScoredDiseases, T: Disease, },
  { ScoredT: ScoredDrugs, T: Drug, },
  { ScoredT: ScoredGenes, T: Gene, },
  { ScoredT: ScoredPathways, T: Pathway, },
  { ScoredT: ScoredPhenotypes, T: Phenotype, },
  { ScoredT: ScoredTissues, T: Tissue, },
].flatMap(({ ScoredT, T }) => [
    MetaNode(`MeanScoredT[${T.label}]`)
      .meta({
        label: 'Mean Scored Values',
        description: 'Take the mean value across multiple scores',
      })
      .inputs({ scored: [ScoredT] })
      .output(ScoredT)
      .resolve(async (props) => {
        const values: Record<string, number[]> = {}
        props.inputs.scored.forEach(scored =>
          scored.forEach(({ term, zscore }) => {
            if (typeof zscore === 'number') {
              if (!(term in values)) values[term] = []
              values[term].push(zscore)
            }
          }))
        return dict.items(values).map(({ key, value }) => ({ term: key, zscore: math.mean(value) }))
      })
      .story(props => `Mean scores were computed.`)
      .build(),
    MetaNode(`AbsMaxScoredT[${T.label}]`)
      .meta({
        label: 'Absolute Max Scored Values',
        description: 'Take the absolute max value across multiple scores',
      })
      .inputs({ scored: [ScoredT] })
      .output(ScoredT)
      .resolve(async (props) => {
        const values: Record<string, number[]> = {}
        props.inputs.scored.forEach(scored =>
          scored.forEach(({ term, zscore }) => {
          if (typeof zscore === 'number') {
            if (!(term in values)) values[term] = []
            values[term].push(zscore)
          }
        }))
        return dict.items(values).map(({ key, value }) => ({ term: key, zscore: math.absmax(value) }))
      })
      .story(props => `Max scores were computed.`)
      .build(),
])
