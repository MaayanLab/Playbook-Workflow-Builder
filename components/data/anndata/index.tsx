import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL, FileC } from '@/components/core/file'
import python from '@/utils/python'
import { z } from 'zod'
import { datafile_icon, file_transfer_icon, transpose_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { downloadUrl } from '@/utils/download'
import { GeneCountMatrix } from '../gene_count_matrix'
import { MetadataMatrix } from '../metadata_matrix'

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
  ))
  .story(props =>
    `The file was parsed as an anndata matrix.`
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
  ))
  .story(props =>
    `An AnnData file was prepared from the input data and metadata.`
  )
  .build()
