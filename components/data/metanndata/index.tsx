import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL, FileC, FilePrompt } from '@/components/core/file'
import python from '@/utils/python'
import { z } from 'zod'
import { datafile_icon, file_icon, file_transfer_icon, transpose_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { downloadUrl } from '@/utils/download'
import { MetaboliteCountMatrix } from '../metabolite_count_matrix'
import { MetadataMatrix } from '../metadata_matrix'
import SafeRender from '@/utils/saferender'

const Matrix = dynamic(() => import('@/app/components/Matrix'))

export const MetAnnData = MetaNode('MetAnnData')
  .meta({
    label: 'Annotated data for metabolites',
    description: 'A metabolite count matrix paired with sample annotations',
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

export const MetAnnDataFromFile = MetaNode('MetAnnDataFromFile')
  .meta({
    label: 'Resolve an MetAnnData Matrix from a File',
    description: 'Ensure a file contains an MetAnnData matrix',
    icon: [file_transfer_icon],
  })
  .inputs({ file: FileURL })
  .output(MetAnnData)
  .resolve(async (props) => await python(
    'components.data.metanndata.metanndata',
    { kargs: [props.inputs.file] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The file${props.inputs && props.inputs.file.description ? ` containing ${props.inputs.file.description}` : ''} was parsed as an metanndata metabolite matrix.`,
    introduction: `The AnnData Python package handles annotated data matrices in memory and on disk\\ref{doi:10.1101/2021.12.16.473007}.`,
    methods: `Metabolite matrices are stored alongside sample metadata annotations in H5 format using AnnData [\\ref{doi:10.1101/2021.12.16.473007}].`,
    legend: `A table showing the basic structure and shape of the uploaded annotated metabolite matrix.`,
  }))
  .build()

export const MetAnnDataFileUpload = MetaNode('MetAnnDataFileUpload')
  .meta({
    label: 'Upload a MetAnnData Matrix',
    description: 'A file containing a MetAnnData matrix',
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
  .output(MetAnnData)
  .prompt(props => <><FilePrompt {...props} />{props.output ? <SafeRender component={MetAnnData.view} props={props.output} /> : null}</>)
  .resolve(async (props) => await python(
    'components.data.metanndata.metanndata',
    { kargs: [props.data] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `A metanndata metabolite matrix${props.data && props.data.description ? ` containing ${props.data.description}` : ''} was uploaded.`,
    introduction: `The AnnData Python package handles annotated data matrices in memory and on disk\\ref{doi:10.1101/2021.12.16.473007}.`,
    methods: `Metabolite matrices are stored alongside sample metadata annotations in H5 format using AnnData [\\ref{doi:10.1101/2021.12.16.473007}].`,
    legend: `A table showing the basic structure and shape of the uploaded annotated metabolite matrix.`,
  }))
  .build()

export const MetAnnDataFromMetaboliteCountMatrixAndMetadataMatrix = MetaNode('MetAnnDataFromMetaboliteCountMatrixAndMetadataMatrix')
  .meta({
    label: 'Resolve an MetAnnData Matrix',
    description: 'Collapse a Metabolite Count Matrix & MetadataMatrix into a single MetAnnData',
    icon: [file_transfer_icon],
  })
  .inputs({ metabolite_count_matrix: MetaboliteCountMatrix, metadata_matrix: MetadataMatrix })
  .output(MetAnnData)
  .resolve(async (props) => await python(
    'components.data.metanndata.metanndata_from_metabolite_count_matrix_and_metadata_matrix',
    { kargs: [props.inputs.metabolite_count_matrix, props.inputs.metadata_matrix] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `An MetAnnData file was prepared from the input data${props.inputs && props.inputs.metabolite_count_matrix.description ? ` containing ${props.inputs.metabolite_count_matrix.description}` : ''} and metadata${props.inputs && props.inputs.metadata_matrix.description ? ` containing ${props.inputs.metadata_matrix.description}` : ''}.`,
    introduction: `The AnnData Python package handles annotated data matrices in memory and on disk\\ref{doi:10.1101/2021.12.16.473007}.`,
    methods: `The metabolite count matrix and metadata matrix were combined into a joint H5 formatted file using AnnData [\\ref{doi:10.1101/2021.12.16.473007}].`,
    legend: `A table showing the basic structure and shape of the uploaded annotated metabolite matrix.`,
  }))
  .build()
