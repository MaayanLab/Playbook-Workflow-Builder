import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL, FileC, FilePrompt } from '@/components/core/file'
import python from '@/utils/python'
import { z } from 'zod'
import { file_icon, file_transfer_icon, metadata_file_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { downloadUrl } from '@/utils/download'
import { clientLoadExample } from './api/meta.tsv/client'
import SafeRender from '@/utils/saferender'

const Matrix = dynamic(() => import('@/app/components/Matrix'))

export const MetadataMatrix = MetaNode('MetadataMatrix')
  .meta({
    label: 'Class Metadata of a Gene Count Matrix',
    description: 'Class metadata for samples in a gene count matrix',
    icon: [metadata_file_icon],
  })
  .codec(FileC.merge(z.object({
    shape: z.tuple([z.number(), z.number()]),
    columns: z.array(z.string()),
    index: z.array(z.string()),
    values: z.array(z.array(z.string())),
    ellipses: z.tuple([z.union([z.number(), z.null()]), z.union([z.number(), z.null()])]),
  })))
  .view(props => {
    return (
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
    )
  })
  .build()

export const MetadataMatrixFromFile = MetaNode('MetadataMatrixFromFile')
  .meta({
    label: 'Resolve a Metadata Matrix from a File',
    description: 'Ensure a file contains a metadata matrix and load it into a standard format',
    icon: [file_transfer_icon],
  })
  .inputs({ file: FileURL })
  .output(MetadataMatrix)
  .resolve(async (props) => await python(
    'components.data.metadata_matrix.metadata_matrix',
    { kargs: [props.inputs.file] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The file${props.inputs && props.inputs.file.description ? ` containing ${props.inputs.file.description}` : ''} was loaded as a metadata matrix.`,
    introduction: `Metadata for gene count matrices provide more information about the cells from which the samples originated, providing context to signatures. Possible metadata categories include sample title, stage of disease, and cell type.`,
    tableLegend: `A table of metadata uploaded by the user.`,
  }))
  .build()

export const MetadataMatrixFileUpload = MetaNode('MetadataMatrixFileUpload')
  .meta({
    label: 'Upload a Metadata Matrix',
    description: 'A file containing a metadata matrix',
    tags: {
      Type: {
        File: 1,
      },
      Cardinality: {
        Matrix: 1,
      },
    },
    icon: [file_icon],
  })
  .codec(FileC)
  .inputs()
  .output(MetadataMatrix)
  .prompt(props => <><FilePrompt {...props} example={clientLoadExample} />{props.output ? <SafeRender component={MetadataMatrix.view} props={props.output} /> : null}</>)
  .resolve(async (props) => await python(
    'components.data.metadata_matrix.metadata_matrix',
    { kargs: [props.data] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `A metadata matrix${props.data && props.data.description ? ` containing ${props.data.description}` : ''} was uploaded.`,
    introduction: `Metadata for gene count matrices provide more information about the cells from which the samples originated, providing context to signatures. Possible metadata categories include sample title, stage of disease, and cell type.`,
    tableLegend: `A table of metadata uploaded by the user.`,
  }))
  .build()
