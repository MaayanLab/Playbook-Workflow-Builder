import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { z } from 'zod'
import { gmt_icon } from '@/icons'
import { Cell, Column, Table } from '@/app/components/Table'
import { GMT } from '../gene_matrix_transpose'
import { downloadBlob } from '@/utils/download'
import { GeneSet } from '@/components/core/set'
import SafeRender from '@/utils/saferender'
import { SELECTED } from '@blueprintjs/core/lib/esm/common/classes'


export const GeneSetCrossing = MetaNode('GeneSetCrossing')
  .meta({
    label: 'Gene Set Crossing',
    description: 'Crossing of multiple gene sets',
    icon: [gmt_icon, gmt_icon],
  })
  .codec(z.object({
    rank:z.number(),
    term1:z.string(),
    term1Length:z.number(),
    term2:z.string(),
    term2Length:z.number(),
    term3:z.string().optional(),
    term3Length:z.number().optional(),
    term4:z.string().optional(),
    term4Length:z.number().optional(),
    term5:z.string().optional(),
    term5Length:z.number().optional(),
    pvalue:z.number(),
    jaccard:z.number(),
    overlap:z.number(),
    genes:z.string()
  }))
  .view(crossing => {
    return (<>
      <div>
        <p><strong>Rank</strong>: {crossing.rank}</p>
        <p><strong>Term 1</strong>: {crossing.term1} ({crossing.term1Length} genes)</p>
        <p><strong>Term 2</strong>: {crossing.term2} ({crossing.term2Length} genes)</p>
        {crossing.term3 && <p><strong>Term 3</strong>: {crossing.term3} ({crossing.term3Length} genes)</p>}
        {crossing.term4 && <p><strong>Term 4</strong>: {crossing.term4} ({crossing.term4Length} genes)</p>}
        {crossing.term5 && <p><strong>Term 5</strong>: {crossing.term5} ({crossing.term5Length} genes)</p>}
        <p><strong>Intersecting Genes</strong>: {crossing.genes} ({crossing.overlap} genes)</p>
      </div>
    </>)
  }).build()

export const GeneSetCrossings = MetaNode('GeneSetCrossings')
  .meta({
    label: 'Gene Set Crossings',
    description: 'Cross 2+ GMTs to identify interesting and unexpected overlaps',
    icon: [gmt_icon, gmt_icon],
  })
  .codec(z.array(z.object({
    rank:z.number(),
    term1:z.string(),
    term1Length:z.number(),
    term2:z.string(),
    term2Length:z.number(),
    term3:z.string().optional(),
    term3Length:z.number().optional(),
    term4:z.string().optional(),
    term4Length:z.number().optional(),
    term5:z.string().optional(),
    term5Length:z.number().optional(),
    pvalue:z.number(),
    jaccard:z.number(),
    overlap:z.number(),
    genes:z.string()
  })))
  .view(crossings => {
    type Row = (typeof crossings)[number]

    const maxTerms = crossings.reduce((max, c) => {
      if (c.term5 !== undefined) return Math.max(max, 5)
      if (c.term4 !== undefined) return Math.max(max, 4)
      if (c.term3 !== undefined) return Math.max(max, 3)
      return Math.max(max, 2)
    }, 2)
  
    const columns: Array<{ name: string; field: keyof Row }> = [
      { name: 'Rank', field: 'rank' },
      { name: 'Term 1', field: 'term1' },
      { name: 'Term 1 Length', field: 'term1Length' },
      { name: 'Term 2', field: 'term2' },
      { name: 'Term 2 Length', field: 'term2Length' },
  
      ...Array.from({ length: maxTerms - 2 }, (_, i) => {
        const n = i + 3
        return [
          { name: `Term ${n}`, field: `term${n}` as keyof Row },
          { name: `Term ${n} Length`, field: `term${n}Length` as keyof Row },
        ]
      }).flat(),
  
      { name: 'p-value', field: 'pvalue' },
      { name: 'Jaccard', field: 'jaccard' },
      { name: 'Overlap', field: 'overlap' },
      { name: 'Genes', field: 'genes' },
    ]

    const formatCell = (field: keyof Row, value: unknown): string => {
      if (value == null) return ''
      if (field === 'pvalue' && typeof value === 'number')
        return value.toExponential(4)
      if (field === 'jaccard' && typeof value === 'number')
        return value.toFixed(6)
      return String(value)
    }
  
    return (
      <Table
        height={500}
        cellRendererDependencies={[crossings]}
        numRows={crossings.length}
        enableGhostCells
        enableFocusedCell
        downloads={{
          JSON: () =>
            downloadBlob(
              new Blob(
                [JSON.stringify(crossings)],
                { type: 'application/json;charset=utf-8' }
              ),
              'data.json'
            ),
  
          TSV: () =>
            downloadBlob(
              new Blob(
                [
                  crossings
                    .map(row =>
                      columns
                        .map(({ field }) => formatCell(field, row[field]))
                        .join('\t')
                    )
                    .join('\n'),
                ],
                { type: 'text/tab-separated-values;charset=utf-8' }
              ),
              'data.tsv'
            ),
        }}
      >
        {columns.map(({ name, field }) => (
          <Column
            key={String(field)}
            name={name}
            cellRenderer={row => (
              <Cell>{formatCell(field, crossings[row][field])}</Cell>
            )}
          />
        ))}
      </Table>
    )
  })
  .build()

export const GMTCrossing = MetaNode('GMTCrossing')
  .meta({
    label: 'GMT Crossing',
    description: 'Cross 2+ GMTs to identify interesting and unexpected overlaps',
    icon: [gmt_icon, gmt_icon],
  })
  .inputs({
    gmts: [GMT],
  })
  .output(GeneSetCrossings)
  .resolve(async (props) => {
    const gmts = props.inputs.gmts.filter((gmt): gmt is NonNullable<typeof gmt> => gmt != null)
    
    return await python(
      'components.data.gene_set_crossing.cross_gmts',
      { kargs: [], kwargs: { gmts:gmts } },
      message => props.notify({ type: 'info', message }),
    )
  })
  .story(props => ({
    abstract: 'GMTs were crossed to identify the combinations of one gene set from each dataset where all gene sets were highly enriched against the intersection of remaining sets.',
    introduction: 'For each combination of one gene set from each dataset, the hypergeometric p-value was computed for each gene set against the intersection of the remaining sets. The most conservative p-value of each crossing was used to rank all combinations, and the top 5000 results were kept.',
    methods: '',
    tableLegend: `A table showing up to the top 5000 crossings between the input GMTs.`,
  }))
  .build()


export const SelectGeneSetCrossing = MetaNode('SelectGeneSetCrossing')
  .meta({
    label: 'Select Gene Set Crossing',
    description: 'Select a crossing from the top combinations of crossed gene sets from multiple libraries',
    icon: [gmt_icon, gmt_icon],
  })
  .codec(GeneSetCrossing.codec)
  .inputs({ crossings: GeneSetCrossings })
  .output(GeneSetCrossing)
  .prompt(props => {
    const { crossings } = props.inputs
    type Row = (typeof crossings)[number]

    const [selected, setSelected] = React.useState<number | null>(
      props.data
        ? crossings.findIndex(c => c.rank === props.data?.rank)
        : null
    )

    const maxTerms: number = crossings.some(c => c.term5 !== undefined) ? 5
      : crossings.some(c => c.term4 !== undefined) ? 4
      : crossings.some(c => c.term3 !== undefined) ? 3
      : 2

    const columns: Array<{ name: string; field: keyof Row }> = [
      { name: 'Rank', field: 'rank' },
      { name: 'Term 1', field: 'term1' },
      { name: 'Term 1 Length', field: 'term1Length' },
      { name: 'Term 2', field: 'term2' },
      { name: 'Term 2 Length', field: 'term2Length' },
      ...Array.from({ length: maxTerms - 2 }, (_, i) => [
        { name: `Term ${i + 3}`, field: `term${i + 3}` as keyof Row },
        { name: `Term ${i + 3} Length`, field: `term${i + 3}Length` as keyof Row },
      ]).flat(),
      { name: 'p-value', field: 'pvalue' },
      { name: 'Jaccard', field: 'jaccard' },
      { name: 'Overlap', field: 'overlap' },
      { name: 'Genes', field: 'genes' },
    ]

    const formatCell = (field: keyof Row, value: unknown): string => {
      if (value == null) return ''
      if (field === 'pvalue' && typeof value === 'number')
        return value.toExponential(4)
      if (field === 'jaccard' && typeof value === 'number')
        return value.toFixed(6)
      return String(value)
    }

    const columnWidths = [
      40,
      ...columns.map(({ name, field }) => {
        const maxLen = crossings.reduce((max, row) => {
          return Math.max(max, formatCell(field, row[field]).length)
        }, name.length)
        return Math.min(Math.max(maxLen * 8 + 24, 60), 300)
      }),
    ]

    return (
      <div className="space-y-2">
        {React.createElement(
          Table,
          {
            cellRendererDependencies: [crossings],
            numRows: crossings.length,
            enableGhostCells: true,
            enableFocusedCell: true,
            onFocusedCell: cell => setSelected(cell.row),
            columnWidths,
            height: 500,
          },
          <Column
            key="select"
            name=""
            cellRenderer={row => (
              <Cell>
                <input
                  type="radio"
                  checked={selected === row}
                  onChange={() => {}}
                  readOnly
                />
              </Cell>
            )}
          />,
          ...columns.map(({ name, field }) => (
            <Column
              key={String(field)}
              name={name}
              cellRenderer={row => (
                <Cell>{formatCell(field, crossings[row][field])}</Cell>
              )}
            />
          ))
        )}

        <button
          className="bp5-button bp5-large"
          disabled={selected === null}
          onClick={() => {
            if (selected !== null) {
              props.submit(crossings[selected]!)
            }
          }}
        >
          Submit
        </button>

        {props.output
          ? <SafeRender component={GeneSetCrossing.view} props={props.output} />
          : null}
      </div>
    )
  })
  .resolve(async props => props.data)
  .story(props => ({
    abstract: `${props.output ? 'The crossing with rank ' + String(props.output.rank) : 'A crossing'} was selected.`,
  }))
  .build()

export const ExtractCrossingTermGeneSet = MetaNode('ExtractCrossingTermGeneSet')
  .meta({
    label: 'Extract Crossing Term Gene Set',
    description: 'Select the term from a crossing and retrieve the corresponding gene set from the crossed GMT',
    icon: [gmt_icon, gmt_icon],
  })
  .codec(z.number())
  .inputs({
    crossing:GeneSetCrossing,
    gmt:GMT
  })
  .output(GeneSet)
  .prompt(props => {
    const [selected, setSelected] = React.useState<number >(
      props.data ?? 1
    )
    const maxTerms = (Object.keys(props.inputs.crossing).length - 5) / 2
    
    return (<>
      <div className="space-y-2">
      
      <select
        className="bp5-select w-full"
        value={selected}
        onChange={e => setSelected(Number(e.target.value))}
      >
        {Array.from({ length: maxTerms }, (_, i) => i + 1).map(i => (
          <option key={i} value={i}>
            Term {i}
          </option>
        ))}
      </select>
        <button
          className="bp5-button bp5-large"
          disabled={selected === null}
          onClick={() => {
            if (selected !== null) {
              props.submit(selected as number)
            }
          }}
        >
          Submit
        </button>
        {props.output
          ? <SafeRender component={GeneSet.view} props={props.output} />
          : null}
      </div>

    </>)
    
  })
  .resolve(async (props) => {
    const key = `term${props.data}` as keyof typeof props.inputs.crossing
    if (props.data < 1 || props.data > 5) {
      throw new Error("Please select a valid term index.")
    }
    const term = props.inputs.crossing[key] as string
    const item = props.inputs.gmt[term]
    if (!item) {
      throw new Error(`Gene set '${term}' was not found in the supplied GMT.`)
    }
    return {
      set:item.set,
      description:term
    }
  })
  .story(props => ({
    abstract: 'The intersecting gene set was extracted from the crossing.',
  }))
  .build()

export const ExtractCrossingGeneSet = MetaNode('ExtractCrossingGeneSet')
  .meta({
    label: 'Extract Crossing Gene Set',
    description: 'Extract the intersecting gene set from a selected gene set crossing',
    icon: [gmt_icon, gmt_icon],
  })
  .inputs({
    crossing:GeneSetCrossing,
  })
  .output(GeneSet)
  .resolve(async (props) => {
    return {
      set:props.inputs.crossing.genes.split(' '),
      description:`Intersecting gene set from ${props.inputs.crossing.term1} & ${props.inputs.crossing.term2}
      ${props.inputs.crossing.term3 ? ' & ' + props.inputs.crossing.term3 : ''}
      ${props.inputs.crossing.term4 ? ' & ' + props.inputs.crossing.term4 : ''}
      ${props.inputs.crossing.term5 ? ' & ' + props.inputs.crossing.term5 : ''} crossing`
    }
  })
  .story(props => ({
    abstract: 'The intersecting gene set was extracted from the crossing.',
  }))
  .build()