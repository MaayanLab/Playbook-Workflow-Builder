import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL } from '@/components/core/file'
import python from '@/utils/python'
import { z } from 'zod'
import { datafile_icon, file_transfer_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { downloadUrl } from '@/utils/download'

const Matrix = dynamic(() => import('@/app/components/Matrix'))

export const GeneCountMatrix = MetaNode('GeneCountMatrix')
  .meta({
    label: 'Gene Count Matrix',
    description: 'A gene count matrix file',
    icon: [datafile_icon],
  })
  .codec(z.object({
    url: z.string(),
    shape: z.tuple([z.number(), z.number()]),
    columns: z.array(z.string()),
    index: z.array(z.string()),
    values: z.array(z.array(z.union([z.number(), z.literal('nan'), z.literal('inf'), z.literal('-inf')]))),
    ellipses: z.tuple([z.union([z.number(), z.null()]), z.union([z.number(), z.null()])]),
  }))
  .view(props => {
    return (
      <div>
        <Matrix
          index={props.index}
          columns={props.columns}
          values={props.values}
          ellipses={props.ellipses}
          shape={props.shape}
          downloads={{
            'URL': () => downloadUrl(props.url)
          }}
        />
      </div>
    )
  })
  .build()

export const GeneCountMatrixFromFile = MetaNode('GeneCountMatrixFromFile')
  .meta({
    label: 'Resolve A Gene Count Matrix from a File',
    description: 'Ensure a file contains a gene count matrix, load it into a standard format',
    icon: [file_transfer_icon],
  })
  .inputs({ file: FileURL })
  .output(GeneCountMatrix)
  .resolve(async (props) => await python(
    'components.data.gene_count_matrix.gene_count_matrix',
    { kargs: [props.inputs.file.url] },
  ))
  .story(props =>
    `The input gene matrix count was converted into a standard format.`
  )
  .build()

  export const Transpose = MetaNode('Transpose')
  .meta({
    label: 'Transpose',
    description: 'A demonstrative transpose operation',
  })
  .inputs({ file: GeneCountMatrix })
  .output(GeneCountMatrix)
  .resolve(async (props) => await python(
    'components.data.gene_count_matrix.transpose',
    { kargs: [props.inputs.file] },
  ))
  .story(props =>
    `The gene count matrix was then transposed`
  )
  .build()
