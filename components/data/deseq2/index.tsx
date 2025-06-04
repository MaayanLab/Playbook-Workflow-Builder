import React from 'react'
import { MetaNode} from '@/spec/metanode';
import python from '@/utils/python'
import { GeneSignature } from '@/components/data/gene_signature'
import { AnnData } from '@/components/data/anndata'

export const DESeq2 = MetaNode('DESeq2')
  .meta({
    label: 'DESeq2 Analysis from AnnData File',
    description: `Construct a table which contains the results of differential gene expression (DGE) analysis between two groups of samples in a dataset. Every row of the table represents a gene; the columns display the results of the differential gene expression analysis.`
  })
  .inputs({
    anndata: AnnData,
  })
  .output(GeneSignature)
  .resolve(async (props) => await python(
    'components.data.deseq2.deseq2_from_matrix',
    { kargs: [props.inputs.anndata] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The AnnData file was then analyzed using differential expression by DESeq2\\ref{doi:10.1186/s13059-014-0550-8}.`,
    introduction: `The AnnData Python package handles annotated data matrices in memory and on disk and can be used to store gene count information. DESeq2 is a pair of normal linear modeilng strategies that is used to analyze read counts from RNA-seq experiments\\ref{doi:10.1186/s13059-014-0550-8}\\ref{doi:10.1093/bioinformatics/btad547}.`,
    methods: `The AnnData package was used to store gene signature data. Differential gene expression analysis was performed using DESeq2\\ref{doi:10.1186/s13059-014-0550-8}\\ref{doi:10.1093/bioinformatics/btad547}.`,
    tableLegend: `A table of results from differential gene expression analysis between two gene sets.`,
  }))
  .build()
