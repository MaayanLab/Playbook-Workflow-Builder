import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL } from '@/components/file'
import python from '@/utils/python'

export const GeneCountMatrix = MetaNode.createData('GeneCountMatrix')
  .meta({
    label: 'Gene Count Matrix',
    description: 'A gene count matrix file',
  })
  .codec<{
    url: string,
    shape: [number, number],
    head_columns: string[], tail_columns: string[],
    head_index: string[], tail_index: string[],
    head_values: number[], tail_values: number[],
  }>()
  .view(props => (
    <div>
      <h2>Gene Count Matrix: {props.url}</h2>
      <span>Shape: ({props.shape[0]}, {props.shape[1]})</span>
      <table>
        <tr>
          <th>&nbsp;</th>
          {props.head_columns.map(col => <th key={col}>{col}</th>)}
          {props.tail_columns.map(col => <th key={col}>{col}</th>)}
        </tr>
        {props.head_index.map((index, i) =>
          <tr key={index}>
            <th>{index}</th>
            {props.head_columns.map((col, j) => <td>{props.head_values[i][j]}</td>)}
            {props.tail_columns.map((col, j) => <td>{props.head_values[i][j]}</td>)}
          </tr>
        )}
        {props.tail_index.map((index, i) =>
          <tr key={index}>
            <th>{index}</th>
            {props.head_columns.map((col, j) => <td>{props.tail_values[i][j]}</td>)}
            {props.tail_columns.map((col, j) => <td>{props.tail_values[i][j]}</td>)}
          </tr>
        )}
      </table>
    </div>
  ))
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
    '@/components/gene-count-matrix/gene_count_matrix.py',
    'gene_count_matrix',
    { kargs: [props.inputs.file] },
  ))
  .build()
