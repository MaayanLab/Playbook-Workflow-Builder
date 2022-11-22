import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL } from '@/components/file'
import python from '@/utils/python'
import { z } from 'zod'
import codecFrom from '@/utils/zod-codec'

export const GeneCountMatrix = MetaNode.createData('GeneCountMatrix')
  .meta({
    label: 'Gene Count Matrix',
    description: 'A gene count matrix file',
  })
  .codec(codecFrom(z.object({
    url: z.string(),
    shape: z.tuple([z.number(), z.number()]),
    columns: z.array(z.string()),
    index: z.array(z.string()),
    values: z.array(z.array(z.union([z.number(), z.literal('nan'), z.literal('inf'), z.literal('-inf')]))),
    ellipses: z.tuple([z.union([z.number(), z.null()]), z.union([z.number(), z.null()])]),
  })))
  .view(props => {
    return (
      <div>
        <h2>Gene Count Matrix: {props.url}</h2>
        <span>Shape: ({props.shape[0]}, {props.shape[1]})</span>
        <table>
          <thead>
            <tr>
              <th>&nbsp;</th>
              {props.columns.flatMap((col, j) => [
                ...(j === props.ellipses[1] ?
                  [<th key={`${j}-ellipse`}>...</th>]
                  : []),
                <th key={j}>{col}</th>,
              ])}
            </tr>
          </thead>
          <tbody>
            {props.index.flatMap((index, i) => [
              ...(i === props.ellipses[0] ?
                [<th key={`${i}-ellipse`}>...</th>]
                : []),
              <tr key={index}>
                <th>{index}</th>
                {props.columns.flatMap((col, j) => [
                  ...(j === props.ellipses[1] ?
                    [<td key={`${j}-ellipse`}>...</td>]
                    : []),
                  <td key={j}>{props.values[i][j]}</td>,
                ])}
              </tr>,
            ])}
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
