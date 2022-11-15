import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL } from '@/components/file'
import python from '@/utils/python'
import * as t from 'io-ts'
import codecFrom from '@/utils/io-ts-codec'

export const GeneCountMatrix = MetaNode.createData('GeneCountMatrix')
  .meta({
    label: 'Gene Count Matrix',
    description: 'A gene count matrix file',
  })
  .codec(codecFrom(t.type({
    url: t.string,
    shape: t.tuple([t.number, t.number]),
    head_columns: t.array(t.string),
    tail_columns: t.array(t.string),
    head_index: t.array(t.string),
    tail_index: t.array(t.string),
    head_values: t.array(t.array(t.union([t.number, t.literal('nan'), t.literal('inf'), t.literal('-inf')]))),
    tail_values: t.array(t.array(t.union([t.number, t.literal('nan'), t.literal('inf'), t.literal('-inf')]))),
  })))
  .view(props => {
    const column_elipse = props.shape[1] > props.head_columns.length + props.tail_columns.length ? '...' : ''
    const index_elipse = props.shape[0] > props.head_index.length + props.tail_index.length ? '...' : ''
    return (
      <div>
        <h2>Gene Count Matrix: {props.url}</h2>
        <span>Shape: ({props.shape[0]}, {props.shape[1]})</span>
        <table>
          <thead>
            <tr>
              <th>&nbsp;</th>
              {props.head_columns.map((col, j) => <th key={j}>{col}</th>)}
              <th>{column_elipse}</th>
              {props.tail_columns.map((col, j) => <th key={j}>{col}</th>)}
            </tr>
          </thead>
          <tbody>
            {props.head_index.map((index, i) =>
              <tr key={index}>
                <th>{index}</th>
                {props.head_columns.map((col, j) => <td key={j}>{props.head_values[i][j]}</td>)}
                <th>{column_elipse}</th>
                {props.tail_columns.map((col, j) => <td key={j}>{props.head_values[i][j]}</td>)}
              </tr>
            )}
            {index_elipse ? (
              <tr>
                <th>...</th>
                {props.head_columns.map((col, j) => <td key={j}>...</td>)}
                <th>{column_elipse}</th>
                {props.tail_columns.map((col, j) => <td key={j}>...</td>)}
              </tr>
            ) : null}
            {props.tail_index.map((index, i) =>
              <tr key={index}>
                <th>{index}</th>
                {props.head_columns.map((col, j) => <td key={j}>{props.tail_values[i][j]}</td>)}
                <th>{column_elipse}</th>
                {props.tail_columns.map((col, j) => <td key={j}>{props.tail_values[i][j]}</td>)}
              </tr>
            )}
          </tbody>
        </table>
      </div>
    )
  })
  .build()

export const GeneCountMatrixFromFile = MetaNode.createProcess('GeneCountMatrixFromFile')
  .meta({
    label: 'Resolve A Gene Count Matrix from a File',
    description: 'Ensure a file contains a gene count matrix, load it into a standard format',
  })
  .codec()
  .inputs({ file: FileURL })
  .output(GeneCountMatrix)
  .resolve(async (props) => await python(
    'components.gene_count_matrix.gene_count_matrix',
    { kargs: [props.inputs.file] },
  ))
  .build()

  export const Transpose = MetaNode.createProcess('Transpose')
  .meta({
    label: 'Transpose',
    description: 'A demonstrative transpose operation',
  })
  .codec()
  .inputs({ file: GeneCountMatrix })
  .output(GeneCountMatrix)
  .resolve(async (props) => await python(
    'components.gene_count_matrix.transpose',
    { kargs: [props.inputs.file] },
  ))
  .build()
