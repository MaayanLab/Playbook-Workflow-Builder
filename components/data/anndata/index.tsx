import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL, FileC, FilePrompt } from '@/components/core/file'
import python from '@/utils/python'
import { z } from 'zod'
import { datafile_icon, file_icon, file_transfer_icon, transpose_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { downloadUrl } from '@/utils/download'
import { GeneCountMatrix } from '../gene_count_matrix'
import { MetadataMatrix } from '../metadata_matrix'
import { clientLoadExample } from  '@/components/data/anndata/api/example.h5ad/client'
import SafeRender from '@/utils/saferender'

const Matrix = dynamic(() => import('@/app/components/Matrix'))

export const AnnData = MetaNode('AnnData')
  .meta({
    label: 'Annotated data',
    description: 'A gene count matrix paired with sample annotations',
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

export const AnnDataFromFile = MetaNode('AnnDataFromFile')
  .meta({
    label: 'Resolve an AnnData Matrix from a File',
    description: 'Ensure a file contains an AnnData matrix',
    icon: [file_transfer_icon],
  })
  .inputs({ file: FileURL })
  .output(AnnData)
  .resolve(async (props) => await python(
    'components.data.anndata.anndata',
    { kargs: [props.inputs.file] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props =>
    `The file${props.inputs && props.inputs.file.description ? ` containing ${props.inputs.file.description}` : ''} was parsed as an anndata matrix.`
  )
  .build()

export const AnnDataFileUpload = MetaNode('AnnDataFileUpload')
  .meta({
    label: 'Upload an AnnData Matrix File',
    description: 'A file containing an AnnData matrix',
    tags: {
      Type: {
        File: 1,
        Gene: 1,
      },
      Cardinality: {
        Matrix: 1,
      },
    },
    icon: [file_icon],
  })
  .codec(FileC)
  .inputs()
  .output(AnnData)
  .prompt(props => <><FilePrompt {...props} example={clientLoadExample} />{props.output ? <SafeRender component={AnnData.view} props={props.output} /> : null}</>)
  .resolve(async (props) => await python(
    'components.data.anndata.anndata',
    { kargs: [props.data] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props =>
    `An anndata matrix${props.data && props.data.description ? ` containing ${props.data.description}` : ''} was uploaded.`
  )
  .build()

export const AnnDataFromGeneCountMatrixAndMetadataMatrix = MetaNode('AnnDataFromGeneCountMatrixAndMetadataMatrix')
  .meta({
    label: 'Resolve an AnnData Matrix',
    description: 'Collapse a Gene Count Matrix & MetadataMatrix into a single AnnData',
    icon: [file_transfer_icon],
  })
  .inputs({ gene_count_matrix: GeneCountMatrix, metadata_matrix: MetadataMatrix })
  .output(AnnData)
  .resolve(async (props) => await python(
    'components.data.anndata.anndata_from_gene_count_matrix_and_metadata_matrix',
    { kargs: [props.inputs.gene_count_matrix, props.inputs.metadata_matrix] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props =>
    `An AnnData file was prepared from the input data${props.inputs && props.inputs.gene_count_matrix.description ? ` containing ${props.inputs.gene_count_matrix.description}` : ''} and metadata${props.inputs && props.inputs.metadata_matrix.description ? ` containing ${props.inputs.metadata_matrix.description}` : ''}.`
  )
  .build()
