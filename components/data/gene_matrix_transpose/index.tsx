import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { Table, Cell, Column, EditableCell } from '@/app/components/Table'
import { GeneSet } from '@/components/core/input/set'
import { downloadBlob } from '@/utils/download'
import dynamic from 'next/dynamic'

const Bp4Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export const GMT = MetaNode(`GMT`)
  .meta({
    label: `Gene Matrix Transpose`,
    description: 'Terms mapped to genes',
  })
  .codec(z.record(z.string(), z.object({ description: z.string().optional(), set: z.array(z.string()) })))
  .view(gmt => {
    const gmt_items = dict.items(gmt)
    return (
      <Table
        height={500}
        cellRendererDependencies={[gmt_items]}
        numRows={gmt_items.length}
        enableGhostCells
        enableFocusedCell
        downloads={{
          JSON: () => downloadBlob(new Blob([JSON.stringify(gmt)], { type: 'application/json;charset=utf-8' }), 'data.json'),
          GMT: () => downloadBlob(new Blob([gmt_items.map(({ key: term, value: { description, set } }) => [term, description||'', ...set].join('\t')).join('\n')], { type: 'text/tab-separated-values;charset=utf-8' }), 'data.gmt'),
        }}
      >
        <Column
          name="Term"
          cellRenderer={row => <Cell key={row+''}>{gmt_items[row].key.toString()}</Cell>}
        />
        <Column
          name="Descrition"
          cellRenderer={row => <Cell key={row+''}>{gmt_items[row].value.description||''}</Cell>}
        />
        <Column
          name="Geneset"
          cellRenderer={row => <Cell key={row+''}>{gmt_items[row].value.set.join('\t')}</Cell>}
        />
      </Table>
    )
  })
  .build()

export const GMTUnion = MetaNode('GMTUnion')
  .meta({
    label: `Compute Union Geneset`,
    description: 'Find the union set of all genes in the GMT'
  })
  .inputs({ gmt: GMT })
  .output(GeneSet)
  .resolve(async (props) => {
    return array.unique(dict.values(props.inputs.gmt).flatMap(({ set: geneset }) => geneset))
  })
  .build()

export const GMTIntersection = MetaNode('GMTIntersection')
  .meta({
    label: `Compute Intersection Geneset`,
    description: 'Find the intersecting set of all genes in the GMT'
  })
  .inputs({ gmt: GMT })
  .output(GeneSet)
  .resolve(async (props) => {
    return dict.values(props.inputs.gmt).reduce(({ set: A }, { set: B }) => ({ set: array.intersection(A, B) })).set
  })
  .build()

export const GMTConsensus = MetaNode('GMTConsensus')
  .meta({
    label: `Compute Consensus Geneset`,
    description: 'Find genes which appear in more than one set'
  })
  .inputs({ gmt: GMT })
  .output(GeneSet)
  .resolve(async (props) => {
    const gene_counts: Record<string, number> = {}
    dict.values(props.inputs.gmt)
      .forEach(({ set: geneset }) =>
        geneset.forEach(gene =>
          gene_counts[gene] = (gene_counts[gene]||0)+1
        )
      )
    return dict.items(gene_counts)
      .filter(({ value }) => value > 1)
      .map(({ key }) => key as string)
  })
  .build()

export const GenesetsToGMT = MetaNode('GenesetsToGMT')
  .meta({
    label: `Assemble GMT from Genesets`,
    description: 'Group multiple independently generated genesets as a single GMT'
  })
  .inputs({ genesets: [GeneSet] })
  .output(GMT)
  .prompt(props => {
    const [terms, setTerms] = React.useState(dict.init(array.arange(props.inputs.genesets.length).map(key => ({ key, value: '' }))))
    const [descriptions, setDescriptions] = React.useState(dict.init(array.arange(props.inputs.genesets.length).map(key => ({ key, value: '' }))))
    React.useEffect(() => {
      if (props.output) {
        setTerms(dict.init(dict.keys(props.output).map((key, i) => ({ key: i, value: key }))))
        setDescriptions(dict.init(dict.values(props.output).map(({ description }, i) => ({ key: i, value: description||'' }))))
      }
    }, [props.inputs, props.output])
    return (
      <div>
        <Table
          height={500}
          cellRendererDependencies={[props.inputs.genesets, terms, descriptions]}
          numRows={props.inputs.genesets.length}
          enableGhostCells
          enableFocusedCell
        >
          <Column
            name="Term"
            cellRenderer={row => <EditableCell
              key={row+''}
              value={terms[row]}
              placeholder={`Geneset from path ${row}`}
              onChange={value => setTerms(terms => ({ ...terms, [row]: value }))}
            />}
          />
          <Column
            name="Descrition"
            cellRenderer={row => <EditableCell
              key={row+''}
              value={descriptions[row]}
              placeholder={`Some optional description`}
              onChange={value => setDescriptions(descriptions => ({ ...descriptions, [row]: value }))}
            />}
          />
          <Column
            name="Geneset"
            cellRenderer={row => <Cell key={row+''}>{props.inputs.genesets[row].join('\t')}</Cell>}
          />
        </Table>
        <Bp4Button
          large
          type="submit"
          text="Submit"
          rightIcon="bring-data"
          onClick={() => props.submit(
            dict.init(
              array.arange(props.inputs.genesets.length)
                .map(i => ({
                  key: terms[i],
                  value: {
                    description: descriptions[i],
                    set: props.inputs.genesets[i],
                  }
                }))
            ))}
        />
      </div>
    )
  })
  .build()

export const GMTConcatenate = MetaNode('GMTConcatenate')
  .meta({
    label: `Concatenate GMTs`,
    description: 'Join several GMTs into one'
  })
  .inputs({ gmts: [GMT] })
  .output(GMT)
  .resolve(async (props) => {
    return dict.init(props.inputs.gmts.flatMap(dict.items))
  })
  .build()
