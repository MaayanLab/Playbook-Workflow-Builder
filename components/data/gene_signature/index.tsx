import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL, FileC, FilePrompt } from '@/components/core/file'
import python from '@/utils/python'
import { z } from 'zod'
import Matrix from '@/app/components/Matrix'
import { datafile_icon, filter_icon, differential_expression_icon, file_transfer_icon, file_icon } from '@/icons'
import { downloadUrl } from '@/utils/download'
import { GMT } from '../gene_matrix_transpose'
import { GeneSet } from '@/components/core/set'
import { ScoredGenes } from '@/components/core/scored'
import SafeRender from '@/utils/saferender'
import { clientLoadExample } from  '@/components/data/gene_signature/api/GTEx_aging_signature_limma.tsv/client'

export const GeneSignature = MetaNode('GeneSignature')
  .meta({
    label: 'Gene Signature',
    description: 'A gene expression signature',
    icon: [differential_expression_icon],
  })
  .codec(FileC.merge(z.object({
    shape: z.tuple([z.number(), z.number()]),
    columns: z.array(z.string()).refine(
      value => value.includes('Pval') && value.includes('AdjPval') && value.includes('LogFC'),
      'Expect columns "Pval" & "AdjPval" & "LogFC"'
    ),
    index: z.array(z.string()),
    values: z.array(z.array(z.union([z.number(), z.literal('inf'), z.literal('-inf')]))),
    ellipses: z.tuple([z.union([z.number(), z.null()]), z.union([z.number(), z.null()])]),
  })))
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
            'URL': () => downloadUrl(props.url, props.filename)
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
    { kargs: [props.inputs.file] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => `The file${props.inputs && props.inputs.file.description ? ` containing ${props.inputs.file.description}` : ''} was loaded as a gene signature.`)
  .build()

export const GeneSigFileUpload = MetaNode('GeneSigFileUpload')
  .meta({
    label: 'Upload a Gene Signature',
    description: 'A table of genes and their significance',
    tags: {
      Type: {
        File: 1,
        Gene: 1,
      },
      Cardinality: {
        Signature: 1,
      },
    },
    icon: [file_icon],
  })
  .codec(FileC)
  .inputs()
  .output(GeneSignature)
  .prompt(props => <><FilePrompt {...props} example={clientLoadExample} />{props.output ? <SafeRender component={GeneSignature.view} props={props.output} /> : null}</>)
  .resolve(async (props) => await python(
    'components.data.gene_signature.gene_signature',
    { kargs: [props.data] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => `A gene signature${props.data && props.data.description ? ` containing ${props.data.description}` : ''} was uploaded.`)
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
    { kargs: [props.inputs.sig] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => `The ${props.inputs && props.inputs.sig.description ? props.inputs.sig.description : 'gene signature'} was reformatted into gene matrix transpose format.`)
  .build()

export const UpGeneSetFromSignature = MetaNode('UpGeneSetFromSignature')
  .meta({
    label: 'Up Gene Set from Signature',
    description: 'Extract significant up-regulated genes from a signature',
    icon: [filter_icon]
  })
  .inputs({ sig: GeneSignature })
  .output(GeneSet)
  .resolve(async (props) => await python (
    'components.data.gene_signature.geneset_from_sig',
    { kargs: [props.inputs.sig, 'up'] }
  ))
  .story(props => `The up-regulated genes were extracted from the ${props.inputs && props.inputs.sig.description ? props.inputs.sig.description : 'gene signature'}.`)
  .build()

export const DownGeneSetFromSignature = MetaNode('DownGeneSetFromSignature')
  .meta({
    label: 'Down Gene Set from Signature',
    description: 'Extract significant down-regulated genes from a signature',
    icon: [filter_icon]
  })
  .inputs({ sig: GeneSignature })
  .output(GeneSet)
  .resolve(async (props) => await python (
    'components.data.gene_signature.geneset_from_sig',
    { kargs: [props.inputs.sig, 'down'] }
  ))
  .story(props => `The down-regulated genes were extracted from the ${props.inputs && props.inputs.sig.description ? props.inputs.sig.description : 'gene signature'}.`)
  .build()

export const ScoredGenesFromSignature = MetaNode('ScoredGenesFromSignature')
  .meta({
    label: 'Scored Genes from Signature',
    description: 'Treat signature as a weighted set of genes',
  })
  .inputs({ sig: GeneSignature })
  .output(ScoredGenes)
  .resolve(async (props) => await python(
    'components.data.gene_signature.scored_genes_from_sig',
    { kargs: [props.inputs.sig] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => `Significant genes were extracted from the ${props.inputs && props.inputs.sig.description ? props.inputs.sig.description : 'gene signature'}.`)
  .build()
