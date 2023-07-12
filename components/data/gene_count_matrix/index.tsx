import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL, FileC } from '@/components/core/file'
import python from '@/utils/python'
import { z } from 'zod'
import { datafile_icon, file_transfer_icon, transpose_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { downloadUrl } from '@/utils/download'

const Matrix = dynamic(() => import('@/app/components/Matrix'))

export const GeneCountMatrix = MetaNode('GeneCountMatrix')
  .meta({
    label: 'Gene Count Matrix',
    description: 'A gene count matrix file',
    icon: [datafile_icon],
  })
  .codec(FileC.merge(z.object({
    shape: z.tuple([z.number(), z.number()]),
    columns: z.array(z.string()),
    index: z.array(z.string()),
    values: z.array(z.array(z.union([z.number(), z.literal('nan'), z.literal('inf'), z.literal('-inf')]))),
    ellipses: z.tuple([z.union([z.number(), z.null()]), z.union([z.number(), z.null()])]),
  })))
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
    label: 'Resolve a Gene Count Matrix from a File',
    description: 'Ensure a file contains a gene count matrix, load it into a standard format',
    icon: [file_transfer_icon],
  })
  .inputs({ file: FileURL })
  .output(GeneCountMatrix)
  .resolve(async (props) => await python(
    'components.data.gene_count_matrix.gene_count_matrix',
    { kargs: [props.inputs.file] },
  ))
  .story(props =>
    `The file${props.inputs && props.inputs.file.description ? ` containing ${props.inputs.file.description}` : ''} was parsed as a gene count matrix.`
  )
  .build()

export const Transpose = MetaNode('Transpose')
  .meta({
    label: 'Transpose',
    description: 'Re-orient the matrix',
    icon: [transpose_icon],
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
