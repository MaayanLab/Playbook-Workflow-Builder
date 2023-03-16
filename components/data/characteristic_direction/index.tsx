import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { MetadataMatrix } from '@/components/data/metadata_matrix'
import { GeneSignature } from '@/components/data/gene_signature'

export const CDSignatureFromCounts = MetaNode('CDSignatureFromCounts')
  .meta({
    label: 'Compute a Characteristic Direction Signature from a Gene Count Matrix',
    description: `From a given counts matrix and metadata, return a charateristic
                  direction differential gene expression signature`
  })
  .inputs({ 
    data: GeneCountMatrix, 
    metadata: MetadataMatrix
  })
  .output(GeneSignature)
  .resolve(async (props) => await python(
    'components.data.characteristic_direction.cd_from_matrix',
    { kargs: [props.inputs.data, props.inputs.metadata]}
  ))
  .build()