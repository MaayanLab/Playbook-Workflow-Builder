import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL } from '@/components/core/file'
import python from '@/utils/python'
import { z } from 'zod'
import { file_transfer_icon, metadata_file_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { downloadUrl } from '@/utils/download'

const Matrix = dynamic(() => import('@/app/components/Matrix'))

export const MetadataMatrix = MetaNode('MetadataMatrix')
  .meta({
    label: 'Class Metadata of a Gene Count Matrix',
    description: 'Class metadata for samples in a gene count matrix',
    icon: [metadata_file_icon],
  })
  .codec(z.object({
    url: z.string(),
    shape: z.tuple([z.number(), z.number()]),
    columns: z.array(z.string()),
    index: z.array(z.string()),
    values: z.array(z.array(z.string())),
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

export const MetadataMatrixFromFile = MetaNode('MetadataMatrixFromFile')
  .meta({
    label: 'Resolve A Metadata Matrix from a File',
    description: 'Ensure a file contains a metadata matrix and load it into a standard format',
    icon: [file_transfer_icon],
  })
  .inputs({ file: FileURL })
  .output(MetadataMatrix)
  .resolve(async (props) => await python(
    'components.data.metadata_matrix.metadata_matrix',
    { kargs: [props.inputs.file] },
  ))
  .build()

// const SampleGroupCol = () => MetaNode('')
//   .meta({
//     label: 'Sample Group Column'
//   })
