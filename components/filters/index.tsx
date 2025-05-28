import React from 'react'
import { z } from 'zod'
import { MetaNode } from '@/spec/metanode'
import { DrugSet, GeneSet, GlycanSet, ProteinSet, VariantSet, RegulatoryElementSet, MetaboliteSet } from '@/components/core/set'
import { DrugRanked, GeneRanked, GlycanRanked, ProteinRanked, VariantRanked, RegulatoryElementRanked, DiseaseRanked, PathwayRanked, PhenotypeRanked, TissueRanked } from '@/components/core/ranked'
import { Disease, Drug, Gene, Glycan, Pathway, Phenotype, Protein, Tissue, Variant, RegulatoryElement, Metabolite } from '@/components/core/primitives'
import { DiseaseTerm, DrugTerm, GeneTerm, GlycanTerm, PathwayTerm, PhenotypeTerm, ProteinTerm, TissueTerm, VariantTerm, RegulatoryElementTerm, MetaboliteTerm } from '@/components/core/term'
import { ScoredDrugs, ScoredGenes, ScoredGlycans, ScoredDiseases, ScoredPathways, ScoredPhenotypes, ScoredProteins, ScoredTissues, ScoredVariants, ScoredRegulatoryElement, ScoredMetabolites } from '@/components/core/scored'
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
  { ScoredT: ScoredVariants, TermT: VariantTerm, T: Variant },
  { ScoredT: ScoredRegulatoryElement, TermT: RegulatoryElementTerm, T: RegulatoryElement },
  { ScoredT: ScoredDiseases, TermT: DiseaseTerm, T: Disease },
  { ScoredT: ScoredPathways, TermT: PathwayTerm, T: Pathway },
  { ScoredT: ScoredPhenotypes, TermT: PhenotypeTerm, T: Phenotype },
  { ScoredT: ScoredTissues, TermT: TissueTerm, T: Tissue },
  { ScoredT: ScoredMetabolites, TermT: MetaboliteTerm, T: Metabolite },
].flatMap(({ ScoredT, TermT, T }) => [
  MetaNode(`TopKScoredT[${ScoredT.spec}]`)
    .meta({
      label: `Top ${ScoredT.meta.label}`,
      description: `Select the top ${ScoredT.meta.label}`,
    })
    .codec(z.object({
      k: z.number(),
    }))
    .inputs({ scored: ScoredT })
    .output(ScoredT)
    .prompt(props => {
      React.useEffect(() => {
        if (!props.data) props.submit({ k: 10 })
      }, [props.data])
      return (
        <div className="flex flex-col gap-2">
          <div className="flex flex-row gap-2">
            {[10, 20, 50, 100].map(k =>
              <button
                key={k}
                className={classNames('btn', {'btn-primary': props.data?.k === k, 'btn-secondary': props.data?.k !== k})}
                onClick={() => {props.submit({ k })}}
              >Top {k}</button>
            )}
          </div>
          {ScoredT.view(props.inputs.scored)}
        </div>
      )
    })
    .resolve(async (props) => {
      return props.inputs.scored.slice(0, props.data.k)
    })
    .story(props => ({
      abstract: `The top ${props.data?.k || 'K'} ${ScoredT.meta.label} were selected.`,
      methods: `The top ${props.data?.k || 'K'} ${pluralize(T.label.toLowerCase())} were selected from ${props.input_refs?.scored}.`,
      legend: `The top ${props.data?.k || 'K'} ${pluralize(T.label.toLowerCase())} from ${props.input_refs?.scored}.`,
    }))
    .build(),
  MetaNode(`OneScoredT[${ScoredT.spec}]`)
    .meta({
      label: `Select One ${TermT.meta.label}`,
      description: `Select one ${TermT.meta.label}`,
    })
    .codec(z.string())
    .inputs({ scored: ScoredT })
    .output(TermT)
    .prompt(props => {
      const scored = props.inputs.scored
      const [selected, setSelected] = React.useState<string | undefined>(props.data)
      React.useEffect(() => {
        if (props.output !== undefined) {
          setSelected(props.output)
        }
      }, [props.output, scored])
      return (
        <Table
          height={500}
          cellRendererDependencies={[scored]}
          rowHeaderCellRenderer={(row) =>
            <div
              className="text-center block"
              onClick={evt => {
                setSelected(scored[row].term)
                props.submit(scored[row].term)
              }}
            >
              <input type="radio" checked={selected === scored[row].term} />
            </div>
          }
          numRows={scored.length}
          enableGhostCells
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
      )
    })
    .resolve(async (props) => {
      if (!props.inputs.scored.some(({ term }) => term === props.data)) {
        throw new Error(`Please select a ${T.label.toLowerCase()} from the table`)
      }
      return props.data
    })
    .story(props => ({
      abstract: props.data ? `${props.data} was chosen for further investigation.` : undefined,
      methods: `${props.data ? props.data : `The ${T.label.toLowerCase()}`} was taken from ${props.input_refs?.scored}.`,
      legend: `${props.data ? props.data : `The ${T.label.toLowerCase()}`} taken from ${props.input_refs?.scored}.`,
    }))
    .build()
])

export const SelectScoredT = [
  { ScoredT: ScoredDiseases },
  { ScoredT: ScoredDrugs },
  { ScoredT: ScoredGenes },
  { ScoredT: ScoredVariants },
  { ScoredT: ScoredRegulatoryElement },
  { ScoredT: ScoredPathways },
  { ScoredT: ScoredPhenotypes },
  { ScoredT: ScoredTissues },
  { ScoredT: ScoredMetabolites },
].flatMap(({ ScoredT }) => [
  MetaNode(`SelectFromScored[${ScoredT.spec}, up]`)
    .meta({
      label: `Up ${ScoredT.meta.label}`,
      description: `Consider only up ${ScoredT.meta.label.toLowerCase()}`,
    })
    .inputs({ scored: ScoredT })
    .output(ScoredT)
    .resolve(async (props) => {
      const scored = props.inputs.scored.filter(({ zscore }) => +zscore > 0)
      scored.sort((a, b) => +b.zscore - +a.zscore)
      return scored
    })
    .story(props => ({
      abstract: `Up ${ScoredT.meta.label.toLowerCase()} were selected.`,
      tableLegend: `A table of ${ScoredT.meta.label.toLowerCase()} with only positive z-scores.`,
    }))
    .build(),
  MetaNode(`SelectFromScored[${ScoredT.spec}, down]`)
    .meta({
      label: `Down ${ScoredT.meta.label}`,
      description: `Consider only down ${ScoredT.meta.label}`,
    })
    .inputs({ scored: ScoredT })
    .output(ScoredT)
    .resolve(async (props) => {
      const scored = props.inputs.scored
        .filter(({ zscore }) => +zscore < 0)
        .map(({ term, zscore }) => ({ term, zscore: -zscore }))
      scored.sort((a, b) => +b.zscore - +a.zscore)
      return scored
    })
    .story(props => ({
      abstract: `Down ${ScoredT.meta.label.toLocaleLowerCase()} were selected.`,
      tableLegend: `A table of ${ScoredT.meta.label.toLowerCase()} with only negative z-scores.`,
    }))
    .build(),
])

export const SetFromScoredT = [
  { ScoredT: ScoredDrugs, SetT: DrugSet, TermT: DrugTerm, T: Drug, },
  { ScoredT: ScoredGenes, SetT: GeneSet, TermT: GeneTerm, T: Gene, },
  { ScoredT: ScoredGlycans, SetT: GlycanSet, TermT: GlycanTerm, T: Glycan, },
  { ScoredT: ScoredProteins, SetT: ProteinSet, TermT: ProteinTerm, T: Protein, },
  { ScoredT: ScoredVariants, SetT: VariantSet, TermT: VariantTerm, T: Variant, },
  { ScoredT: ScoredRegulatoryElement, SetT: RegulatoryElementSet, TermT: RegulatoryElementTerm, T: RegulatoryElement, },
  { ScoredT: ScoredMetabolites, SetT: MetaboliteSet, TermT: MetaboliteTerm, T: Metabolite, },
].flatMap(({ ScoredT, SetT, TermT, T }) => [
  MetaNode(`SetFromScored[${ScoredT.spec}]`)
    .meta({
      label: `Set from ${ScoredT.meta.label}`,
      description: `Construct set with ${ScoredT.meta.label.toLowerCase()}`,
    })
    .inputs({ scored: ScoredT })
    .output(SetT)
    .resolve(async (props) => ({ set: props.inputs.scored.map(({ term }) => term) }))
    .story(props => ({
      abstract: `A set was constructed using the ${ScoredT.meta.label.toLowerCase()}.`,
      tableLegend: `The set of drugs from ${props.input_refs?.scored}.`,
    }))
    .build(),
  MetaNode(`OneSetT[${SetT.spec}]`)
    .meta({
      label: `Select One ${TermT.meta.label}`,
      description: `Select one ${TermT.meta.label}`,
    })
    .codec(z.string())
    .inputs({ set: SetT })
    .output(TermT)
    .prompt(props => {
      const set = props.inputs.set.set
      const [selected, setSelected] = React.useState<string | undefined>(props.data)
      React.useEffect(() => {
        if (props.output !== undefined) {
          setSelected(props.output)
        }
      }, [props.output, set])
      return (
        <Table
          height={500}
          cellRendererDependencies={[set]}
          rowHeaderCellRenderer={(row) =>
            <div
              className="text-center block"
              onClick={evt => {
                setSelected(set[row])
                props.submit(set[row])
              }}
            >
              <input type="radio" checked={selected === set[row]} />
            </div>
          }
          numRows={set.length}
          enableGhostCells
        >
          <Column
            name={T.label}
            cellRenderer={row => <Cell key={row+''}>{set[row]}</Cell>}
          />
        </Table>
      )
    })
    .resolve(async (props) => {
      if (!props.inputs.set.set.some((term) => term === props.data)) {
        throw new Error(`Please select a ${T.label.toLowerCase()} from the table`)
      }
      return props.data
    })
    .story(props => ({
      abstract: props.data ? `${props.data} was chosen for further investigation.` : undefined,
      tableLegend: `${props.data ? props.data : `The ${T.label.toLowerCase()}`} taken from ${props.input_refs?.set}.`,
    }))
    .build(),
  MetaNode(`SomeSetT[${SetT.spec}]`)
    .meta({
      label: `Select Some ${pluralize(TermT.meta.label)}`,
      description: `Select some ${pluralize(TermT.meta.label)}`,
    })
    .codec(z.record(z.string(), z.literal(true)))
    .inputs({ set: SetT })
    .output(SetT)
    .prompt(props => {
      const set = props.inputs.set.set
      const [selected, setSelected] = React.useState(props.data ? props.data : {} as Record<string, true>)
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
        <>
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
          >
            <Column
              name={T.label}
              cellRenderer={row => <Cell key={row+''}>{set[row]}</Cell>}
            />
          </Table>
          <button className="bp5-button bp5-large" onClick={async () => {
            props.submit(selected)
          }}>Submit</button>
        </>
      )
    })
    .resolve(async (props) => {
      const set = new Set()
      props.inputs.set.set.forEach(item => set.add(item))
      if (Object.keys(props.data).some(key => !set.has(key))) {
        throw new Error('Please select genes from the table')
      }
      return {
        description: props.inputs.set.description ? `Filtered ${props.inputs.set.description}` : undefined,
        set: props.inputs.set.set.filter((item) => !props.data[item]),
      }
    })
    .story(props =>  ({
      abstract: `Some ${pluralize(T.label.toLowerCase())} were selected for further investigation.`,
      tableLegend: `Some ${pluralize(T.label.toLowerCase())} taken from ${props.input_refs?.set}.`,
    }))
    .build(),
])

export const ReduceMultiScoredT = [
  { ScoredT: ScoredDiseases, T: Disease, },
  { ScoredT: ScoredDrugs, T: Drug, },
  { ScoredT: ScoredGenes, T: Gene, },
  { ScoredT: ScoredVariants, T: Variant, },
  { ScoredT: ScoredRegulatoryElement, T: RegulatoryElement, },
  { ScoredT: ScoredPathways, T: Pathway, },
  { ScoredT: ScoredPhenotypes, T: Phenotype, },
  { ScoredT: ScoredTissues, T: Tissue, },
  { ScoredT: ScoredMetabolites, T: Metabolite, },
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
        const results = dict.items(values).map(({ key, value }) => ({ term: key, zscore: math.mean(value) }))
        results.sort((a, b) => b.zscore - a.zscore)
        return results
      })
      .story(props => ({
        abstract: `The mean across multiple ${ScoredT.meta.label} is computed.`,
        tableLegend: `The mean across${props.input_refs?.scored ? ` ${(props.input_refs?.scored as string[]).join(', ')}` : ''}.`,
      }))
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
        const results = dict.items(values).map(({ key, value }) => ({ term: key, zscore: math.absmax(value) }))
        results.sort((a, b) => b.zscore - a.zscore)
        return results
      })
      .story(props => ({
        abstract: `The absolute maximum across multiple ${ScoredT.meta.label} is computed.`,
        tableLegend: `The absolute maximum${props.input_refs?.scored ? ` across ${(props.input_refs.scored as string[]).join(', ')}` : ''}.`,
      }))
      .build(),
])


export const TopKRankedT = [
  { RankedT: DrugRanked, TermT: DrugTerm, T: Drug },
  { RankedT: GeneRanked, TermT: GeneTerm, T: Gene },
  { RankedT: VariantRanked, TermT: VariantTerm, T: Variant },
  { RankedT: RegulatoryElementRanked, TermT: RegulatoryElementTerm, T: RegulatoryElement },
  { RankedT: DiseaseRanked, TermT: DiseaseTerm, T: Disease },
  { RankedT: PathwayRanked, TermT: PathwayTerm, T: Pathway },
  { RankedT: PhenotypeRanked, TermT: PhenotypeTerm, T: Phenotype },
  { RankedT: TissueRanked, TermT: TissueTerm, T: Tissue },
].flatMap(({ RankedT, TermT, T }) => [
  MetaNode(`TopKRankedT[${RankedT.spec}]`)
    .meta({
      label: `Top ${RankedT.meta.label}`,
      description: `Select the top ${RankedT.meta.label}`,
    })
    .codec(z.object({
      k: z.number(),
    }))
    .inputs({ ranked: RankedT })
    .output(RankedT)
    .prompt(props => {
      React.useEffect(() => {
        if (!props.data) props.submit({ k: 10 })
      }, [props.data])
      return (
        <div className="flex flex-col gap-2">
          <div className="flex flex-row gap-2">
            {[10, 20, 50, 100].map(k =>
              <button
                key={k}
                className={classNames('btn', {'btn-primary': props.data?.k === k, 'btn-secondary': props.data?.k !== k})}
                onClick={() => {props.submit({ k })}}
              >Top {k}</button>
            )}
          </div>
          {RankedT.view(props.inputs.ranked)}
        </div>
      )
    })
    .resolve(async (props) => {
      return {
        ranked: props.inputs.ranked.ranked.slice(0, props.data.k),
        description: `Top ${props.data.k} ${props.inputs.ranked.description}`,
      }
    })
    .story(props => ({
      abstract: `The top ${props.data?.k || 'K'} ${RankedT.meta.label} were selected.`,
      tableLegend: `The top ${props.data?.k || 'K'} drugs from ${props.input_refs?.ranked}.`,
    }))
    .build(),
  MetaNode(`OneRankedT[${RankedT.spec}]`)
    .meta({
      label: `Select One ${TermT.meta.label}`,
      description: `Select one ${TermT.meta.label}`,
    })
    .codec(z.string())
    .inputs({ ranked: RankedT })
    .output(TermT)
    .prompt(props => {
      const ranked = props.inputs.ranked.ranked
      const [selected, setSelected] = React.useState<string | undefined>(props.data)
      React.useEffect(() => {
        if (props.output !== undefined) {
          setSelected(props.output)
        }
      }, [props.output, ranked])
      return (
        <>
          <Table
            height={500}
            cellRendererDependencies={[ranked]}
            rowHeaderCellRenderer={(row) =>
              <div
                className="text-center block"
                onClick={evt => {
                  setSelected(ranked[row])
                  props.submit(ranked[row])
                }}
              >
                <input type="radio" checked={selected === ranked[row]} />
              </div>
            }
            numRows={ranked.length}
            enableGhostCells
          >
            <Column
              name={T.label}
              cellRenderer={row => <Cell key={row+''}>{ranked[row]}</Cell>}
            />
          </Table>
        </>
      )
    })
    .resolve(async (props) => {
      if (!props.inputs.ranked.ranked.some((term) => term === props.data)) {
        throw new Error(`Please select a ${T.label.toLowerCase()} from the table`)
      }
      return props.data
    })
    .story(props => ({
      abstract: props.data ? `${props.data} was chosen for further investigation.` : undefined,
      legend: `${props.data ? props.data : `The ${T.label.toLowerCase()}`} taken from ${props.input_refs?.ranked}.`,
    }))
    .build()
])
