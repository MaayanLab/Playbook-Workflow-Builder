import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL, FileC, FilePrompt } from '@/components/core/file'
import python from '@/utils/python'
import { z } from 'zod'
import { datafile_icon, file_icon, file_transfer_icon, transpose_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { downloadUrl } from '@/utils/download'
import { clientLoadExample } from  '@/components/data/gene_count_matrix/api/example.tsv/client'
import SafeRender from '@/utils/saferender'

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
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The file${props.inputs && props.inputs.file.description ? ` containing ${props.inputs.file.description}` : ''} was parsed as a gene count matrix.`,
    tableLegend: `A table showing the basic structure and shape of the uploaded gene count matrix. Rows represent genes, columns represent samples, and values show the number of mapped reads.`,
  }))
  .build()

export const GeneCountMatrixFileUpload = MetaNode('GeneCountMatrixFileUpload')
  .meta({
    label: 'Upload a Gene Count Matrix',
    description: 'A file containing a gene count matrix',
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
  .output(GeneCountMatrix)
  .prompt(props => <><FilePrompt {...props} example={clientLoadExample} />{props.output ? <SafeRender component={GeneCountMatrix.view} props={props.output} /> : null}</>)
  .resolve(async (props) => await python(
    'components.data.gene_count_matrix.gene_count_matrix',
    { kargs: [props.data] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The gene count matrix${props.data && props.data.description ? ` containing ${props.data.description}` : ''} was uploaded.`,
    tableLegend: `A table showing the basic structure and shape of the uploaded gene count matrix. Rows represent genes, columns represent samples, and values show the number of mapped reads.`,
  }))
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
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The gene count matrix was then transposed.`,
    methods: `The gene count matrix in ${props.input_refs?.file} as transposed to produce ${props.output_ref}.`,
    tableLegend: `A table showing the basic structure and shape of the transposed gene count matrix. Rows represent columns, columns represent genes, and values show the number of mapped reads.`,
  }))
  .build()
