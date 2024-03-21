import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL, FileC, FilePrompt } from '@/components/core/file'
import python from '@/utils/python'
import { z } from 'zod'
import { datafile_icon, file_icon, file_transfer_icon, transpose_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { downloadUrl } from '@/utils/download'

const Matrix = dynamic(() => import('@/app/components/Matrix'))

export const MetaboliteCountMatrix = MetaNode('MetaboliteCountMatrix')
.meta({
    label: 'Metabolite Count Matrix',
    description: 'A metabolite count matrix file',
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
            'URL': () => downloadUrl(props.url, props.filename)
          }}
        />
      </div>
    )
})
.build()

export const MetaboliteCountMatrixFromFile = MetaNode('MetaboliteCountMatrixFromFile')
.meta({
    label: 'Resolve a Metabolite Count Matrix from a File',
    description: 'Ensure a file contains a metabolite count matrix, load it into a standard format',
    icon: [file_transfer_icon],
})
.inputs({ file: FileURL })
.output(MetaboliteCountMatrix)
.resolve(async (props) => await python(
  'components.data.metabolite_count_matrix.metabolite_count_matrix',
   { kargs: [props.inputs.file] },
   message => props.notify({ type: 'info', message }),
))
.story(props =>
  `The file${props.inputs && props.inputs.file.description ? ` containing ${props.inputs.file.description}` : ''} was parsed as a gene count matrix.`
)
.build()


export const MetaboliteCountMatrixFileUpload = MetaNode('MetaboliteCountMatrixFileUpload')
  .meta({
    label: 'Upload a Metabolite Count Matrix',
    description: 'A file containing a metabolite count matrix',
    tags: {
      Type: {
        File: 1,
        Metabolite: 1,
      },
      Cardinality: {
        Matrix: 1,
      },
    },
    icon: [file_icon],
  })
  .codec(FileC)
  .inputs()
  .output(MetaboliteCountMatrix)
  .prompt(props => <FilePrompt {...props} />)
  .resolve(async (props) => await python(
    'components.data.metabolite_count_matrix.metabolite_count_matrix',
    { kargs: [props.data] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props =>
    `A metabolite count matrix${props.data && props.data.description ? ` containing ${props.data.description}` : ''} was uploaded.`
  )
  .build()
