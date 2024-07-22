import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneSignature } from '@/components/data/gene_signature'
import { AnnData } from '@/components/data/anndata'
import Citable from '@/utils/citations'

export const CDSignatureFromCounts = MetaNode('CDSignatureFromCounts')
  .meta({
    label: 'Compute a Characteristic Direction Signature from a Gene Count Matrix',
    description: `From a given counts matrix and metadata, return a Charateristic
                  Direction differential gene expression signature`
  })
  .inputs({
    anndata: AnnData,
  })
  .output(GeneSignature)
  .resolve(async (props) => await python(
    'components.data.characteristic_direction.cd_from_matrix',
    { kargs: [props.inputs.anndata] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({ abstract: Citable.text`Characteristic direction [${Citable.doi('10.1186/1471-2105-15-79')}] is applied to the input anndata${props.inputs && props.inputs.anndata.description ? ` containing ${props.inputs.anndata.description}` : ''}.` }))
  .build()
