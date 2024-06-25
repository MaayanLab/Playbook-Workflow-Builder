import React from 'react'
import { MetaNode} from '@/spec/metanode';
import python from '@/utils/python'
import { GeneSignature } from '@/components/data/gene_signature'
import { AnnData } from '@/components/data/anndata'

export const TTest = MetaNode('TTest')
  .meta({
    label: 'T-Test Differential Expression of AnnData File',
    description: `Construct a table which contains the results of differential gene expression (DGE) analysis between two groups of samples in a dataset. Every row of the table represents a gene; the columns display the results of the differential gene expression analysis.`
  })
  .inputs({
    anndata: AnnData,
  })
  .output(GeneSignature)
  .resolve(async (props) => await python(
    'components.data.ttest.ttest_from_matrix',
    { kargs: [props.inputs.anndata] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The AnnData file was then analyzed using differential expression by Students T-Test.`,
  }))
  .build()
