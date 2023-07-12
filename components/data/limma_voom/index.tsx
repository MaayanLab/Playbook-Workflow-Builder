import React from 'react'
import { MetaNode} from '@/spec/metanode';
import python from '@/utils/python'
import { GeneSignature } from '@/components/data/gene_signature'
import { AnnData } from '@/components/data/anndata'

export const limmavoom = MetaNode('limmavoom')
  .meta({
    label: 'limmavoom',
    description: `limmavoom`
  })
  .inputs({
    anndata: AnnData,
  })
  .output(GeneSignature)
  .resolve(async (props) => await python(
    'components.data.limma_voom.limma_voom_from_matrix',
    { kargs: [props.inputs.anndata]}
  ))
  .story(props =>
    `The limmavoom plug-in embeds a table which contains the results of differential gene expression (DGE) analysis between two groups of samples in a dataset. Every row of the table represents a gene; the columns display the results of the differential gene expression analysis.`
  )
  .build()
 

