import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL } from '@/components/core/file'
import python from '@/utils/python'
import { z } from 'zod'
import Matrix from '@/app/components/Matrix'
import { downloadUrl } from '@/utils/download'
import { GMT } from '../gene_matrix_transpose'

export const GeneSignature = MetaNode('GeneSignature')
  .meta({
    label: 'Gene Signature',
    description: 'A gene expression signature',
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
  })
  .inputs({ file: FileURL })
  .output(GeneSignature)
  .resolve(async (props) => await python(
    'components.data.gene_signature.gene_signature',
    { kargs: [props.inputs.file] },
  ))
  .build()

export const GMTFromSignature = MetaNode('GMTFromSignature')
  .meta({
    label: 'GMT From Signature',
    description: 'Generate a GMT of up and down genes from a signature'
  })
  .inputs({ sig: GeneSignature })
  .output(GMT)
  .resolve(async (props) => await python(
    'components.data.gene_signature.gmt_from_sig',
    { kargs: [props.inputs.sig] } 
  ))
  .build()

