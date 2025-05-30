import React from 'react'
import { MetaNode } from '@/spec/metanode'
import python from '@/utils/python'
import { GeneSignature } from '@/components/data/gene_signature'
import { AnnData } from '@/components/data/anndata'

export const CDSignatureFromCounts = MetaNode('CDSignatureFromCounts')
  .meta({
    label: 'Compute a Characteristic Direction Signature from a Gene Count Matrix',
    description: `From a given counts matrix and metadata, return a Charateristic Direction differential gene expression signature`,
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
  .story(props => ({
    abstract: `Characteristic direction\\ref{doi:10.1186/1471-2105-15-79} is applied to the input anndata${props.inputs && props.inputs.anndata.description ? ` containing ${props.inputs.anndata.description}` : ''}.`,
    introduction: `The AnnData Python package handles annotated data matrices in memory and on disk\\ref{doi:10.1101/2021.12.16.473007}.`,
    methods: `Characteristic direction\\ref{doi:10.1186/1471-2105-15-79} is applied to the input anndata.`,
    tableLegend: `A table of significantly deferentially expressed genes and characteristic direction coefficients.`,
  }))
  .build()
