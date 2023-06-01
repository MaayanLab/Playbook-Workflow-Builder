import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL } from '@/components/core/file'
import python from '@/utils/python'
import { z } from 'zod'
import Matrix from '@/app/components/Matrix'
import { datafile_icon, filter_icon, differential_expression_icon, file_transfer_icon } from '@/icons'
import { downloadUrl } from '@/utils/download'
import { GMT } from '../gene_matrix_transpose'
import { GeneSet } from '@/components/core/input/set'

export const GeneSignature = MetaNode('GeneSignature')
  .meta({
    label: 'Gene Signature',
    description: 'A gene expression signature',
    icon: [differential_expression_icon],
  })
  .codec(z.object({
    url: z.string(),
    shape: z.tuple([z.number(), z.number()]),
    columns: z.array(z.string()),
    index: z.array(z.string()),
    values: z.array(z.array(z.union([z.number(), z.literal('inf'), z.literal('-inf')]))),
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

export const GeneSigFromFile = MetaNode('GeneSigFromFile')
  .meta({
    label: 'Resolve a Gene Signature from a File',
    description: 'Load a gene signature into a standard format',
    icon: [file_transfer_icon]
  })
  .inputs({ file: FileURL })
  .output(GeneSignature)
  .resolve(async (props) => await python(
    'components.data.gene_signature.gene_signature',
    { kargs: [props.inputs.file.url] },
  ))
  .build()

export const GMTFromSignature = MetaNode('GMTFromSignature')
  .meta({
    label: 'GMT from Signature',
    description: 'Generate a GMT of up and down genes from a signature'
  })
  .inputs({ sig: GeneSignature })
  .output(GMT)
  .resolve(async (props) => await python(
    'components.data.gene_signature.gmt_from_sig',
    { kargs: [props.inputs.sig] }
  ))
  .build()

export const UpGeneSetFromSignature = MetaNode('UpGeneSetFromSignature')
  .meta({
    label: 'Up Gene Set from Signature',
    description: 'Extract top 250 up-regulated genes from a signature',
    icon: [filter_icon]
  })
  .inputs({ sig: GeneSignature })
  .output(GeneSet)
  .resolve(async (props) => await python (
    'components.data.gene_signature.geneset_from_sig',
    { kargs: [props.inputs.sig, 'up'] }
  ))
  .build()

export const DownGeneSetFromSignature = MetaNode('DownGeneSetFromSignature')
  .meta({
    label: 'Down Gene Set from Signature',
    description: 'Extract top 250 down-regulated genes from a signature',
    icon: [filter_icon]
  })
  .inputs({ sig: GeneSignature })
  .output(GeneSet)
  .resolve(async (props) => await python (
    'components.data.gene_signature.geneset_from_sig',
    { kargs: [props.inputs.sig, 'down'] }
  ))
  .build()
