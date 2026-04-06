import { DrugSet, GeneSet } from "@/components/core/set"
import { ScoredGenes } from "@/components/core/scored"
import dynamic from 'next/dynamic'
import python from '@/utils/python'
import { downloadBlob } from '@/utils/download'
import { perturbseqr_icon } from "@/icons"
import { MetaNode } from "@/spec/metanode"
import * as array from '@/utils/array'
import { z } from 'zod'

const perturbseqr_url = 'https://perturbseqr.maayanlab.cloud'
const Matrix = dynamic(() => import('@/app/components/Matrix'))

const PerturbSeqrValue = z.union([z.string(), z.number(), z.boolean(), z.null()])

export const PerturbSeqrEnrichmentResults = MetaNode(`Perturb-SeqrEnrichmentResults`)
  .meta({
    label: 'Perturb-Seqr Enrichment Results',
    description: 'Enrichment results from Perturb-Seqr',
    icon: [perturbseqr_icon],
    external: true,
  })
  .codec(z.object({
    shape: z.tuple([z.number(), z.number()]),
    index: z.array(z.string()),
    columns: z.array(z.string()),
    values: z.array(z.array(PerturbSeqrValue)),
    ellipses: z.tuple([z.number().nullable(), z.number().nullable()]),
  }))
  .view(results => (
    <Matrix
      index={results.index}
      columns={results.columns}
      values={results.values.map(row => row.map(v => String(v ?? '')))}
      ellipses={results.ellipses}
      shape={results.shape}
      downloads={{
        JSON: () => downloadBlob(
          new Blob([JSON.stringify(results)], { type: 'application/json;charset=utf-8' }),
          'perturbseqr_enrichment.json'
        ),
        CSV: () => downloadBlob(
          new Blob([
            [
              ['', ...results.columns].join(','),
              ...results.index.map((idx, i) =>
                [idx, ...results.values[i].map(v => v ?? '')].join(',')
              ),
            ].join('\n')
          ], { type: 'text/csv;charset=utf-8' }),
          'perturbseqr_enrichment.csv'
        ),
      }}
    />
  ))
  
  .build()

export const PerturbSeqrGeneSet = MetaNode(`Perturb-SeqrGeneSet`)
  .meta({
    label: 'Perturb-Seqr Gene Set',
    description: 'A gene set uploaded to Perturb-Seqr',
    icon: [perturbseqr_icon],
    external: true,
  })
  .codec(z.object({
    id: z.string(),
  }))
  .view(props => (
    <div className="prose max-w-none">
      <p>Gene set uploaded to Perturb-Seqr with ID: <code>{props.id}</code></p>
      <p>Use <em>Fetch Perturb-Seqr Enrichment Results</em> to load the results table.</p>
    </div>
  ))
  .build()

export const PerturbSeqrGeneSignature = MetaNode(`Perturb-SeqrGeneSignature`)
  .meta({
    label: 'Perturb-Seqr Gene Signature',
    description: 'A gene signature uploaded to Perturb-Seqr',
    icon: [perturbseqr_icon],
    external: true,
  })
  .codec(z.object({ up_id: z.string(), down_id: z.string() }))
  .view(props => (
    <div className="prose max-w-none">
      <p>Gene signature uploaded to Perturb-Seqr.</p>
      <p>Up ID: <code>{props.up_id}</code> / Down ID: <code>{props.down_id}</code></p>
      <p>Use <em>Fetch Perturb-Seqr Enrichment Results</em> to load the results table.</p>
    </div>
  ))
  .build()

export const PerturbSeqrEnrichmentAnalysis = MetaNode(`Perturb-SeqrEnrichmentAnalysis`)
  .meta({
    label: 'LINCS L1000 Signature Search (Perturb-Seqr)',
    description: 'Use Perturb-Seqr to identify small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to the gene set',
    icon: [perturbseqr_icon],
    external: true,
  })
  .inputs({ gene_set: GeneSet })
  .output(PerturbSeqrEnrichmentResults)
  .resolve(async (props) => {
    const req = await fetch(`${perturbseqr_url}/graphql`, {
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
      },
      method: 'POST',
      body: JSON.stringify({"operationName":"AddUserGeneSet","variables":{"description":`${props.inputs.gene_set.description ?? 'Gene set'} from pwb`,"genes":props.inputs.gene_set.set},"query":"mutation AddUserGeneSet($genes: [String], $description: String = \"\") {\n  addUserGeneSet(input: {genes: $genes, description: $description}) {\n    userGeneSet {\n      id\n      __typename\n    }\n    __typename\n  }\n}\n"}),
    })
    if (!req.ok) throw new Error('Failed to submit gene set to Perturb-Seqr')
    const { data: { addUserGeneSet: { userGeneSet: { id } } } } = await req.json()
    return await python(
      'components.service.perturbseqr.fetch_geneset_enrichment_results',
      { kargs: [id] },
      message => props.notify({ type: 'info', message }),
    )
  })
  .story(props => ({
    abstract: `Small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to gene sets with ${props.inputs?.gene_set?.description ? props.inputs.gene_set.description : 'the gene set'} were identified using Perturb-Seqr.`,
    legend: `A listing of small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to gene sets with ${props.inputs?.gene_set?.description ? props.inputs.gene_set.description : 'the gene set'} found using Perturb-Seqr.`,
  }))
  .build()

export const PerturbSeqrSignatureEnrichmentAnalysis = MetaNode(`Perturb-SeqrSignatureEnrichmentAnalysis`)
  .meta({
    label: 'Perturb-Seqr Signature Search',
    description: 'Use Perturb-Seqr to identify small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to the gene signature',
    icon: [perturbseqr_icon],
    external: true,
  })
  .inputs({ scored_genes: ScoredGenes })
  .output(PerturbSeqrEnrichmentResults)
  .resolve(async (props) => {
    const postSet = async (description: string, genes: string[]) => {
      const req = await fetch(`${perturbseqr_url}/graphql`, {
        headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
        method: 'POST',
        body: JSON.stringify({
          operationName: 'AddUserGeneSet',
          variables: { description, genes },
          query: 'mutation AddUserGeneSet($genes: [String], $description: String = "") {\n  addUserGeneSet(input: {genes: $genes, description: $description}) {\n    userGeneSet {\n      id\n      __typename\n    }\n    __typename\n  }\n}\n',
        }),
      })
      if (!req.ok) throw new Error('Failed to submit gene set to Perturb-Seqr')
      const { data: { addUserGeneSet: { userGeneSet: { id } } } } = await req.json()
      return id
    }

    const [up_id, down_id] = await Promise.all([
      postSet('Up gene set from pwb', props.inputs.scored_genes
        .filter(g => g.zscore === 'inf' || (typeof g.zscore === 'number' && g.zscore > 0))
        .map(g => g.term)
        .slice(0,500)),
      postSet('Down gene set from pwb', props.inputs.scored_genes
        .filter(g => g.zscore === '-inf' || (typeof g.zscore === 'number' && g.zscore < 0))
        .map(g => g.term)
        .slice(-500)),
    ])

    return await python(
      'components.service.perturbseqr.fetch_enrichment_results',
      { kargs: [up_id, down_id] },
      message => props.notify({ type: 'info', message }),
    )
  })
  .story(props => ({
    abstract: `Small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to the signature were identified using Perturb-Seqr.`,
    legend: `A table of drug and gene perturbation enrichment results from Perturb-Seqr.`,
  }))
  .build()

function makePerturbSeqrExtractNode<T>(
  perturbationType: 'Gene' | 'Drug',
  direction: 'Mimickers' | 'Reversers',
  SetT: typeof DrugSet | typeof GeneSet
) {
  return MetaNode(`ExtractPerturb-Seqr${perturbationType}${direction}`)
    .meta({
      label: `Extract ${perturbationType} ${direction} from Perturb-Seqr Results`,
      description: `Extract ${perturbationType.toLocaleLowerCase()} ${direction.toLocaleLowerCase()} from Perturb-Seqr enrichment results`,
      icon: [perturbseqr_icon],
      external: true,
    })
    .inputs({ searchResults: PerturbSeqrEnrichmentResults })
    .output(SetT)
    .resolve(async (props) => {
      return await python(
        'components.service.perturbseqr.extract_perturbation_set',
        { kargs: [props.inputs.searchResults, perturbationType, direction] },
        message => props.notify({ type: 'info', message }),
      )
    })
    .story(props => ({
      abstract: `${perturbationType} ${direction.toLocaleLowerCase()} were extracted from the Perturb-Seqr enrichment results.`,
      methods: `Rows whose dataset matched the ${perturbationType.toLowerCase()} perturbation libraries with p-value < 0.05 were selected.`,
      tableLegend: `A set of ${perturbationType.toLocaleLowerCase()} ${direction.toLocaleLowerCase()} identified from the Perturb-Seqr enrichment results.`,
    }))
    .build()
}

export const ExtractPerturbSeqrMimickerGenes = makePerturbSeqrExtractNode('Gene', 'Mimickers', GeneSet)
export const ExtractPerturbSeqrMimickerDrugs = makePerturbSeqrExtractNode('Drug', 'Mimickers', DrugSet)
export const ExtractPerturbSeqrReverserGenes = makePerturbSeqrExtractNode('Gene', 'Reversers', GeneSet)
export const ExtractPerturbSeqrReverserDrugs = makePerturbSeqrExtractNode('Drug', 'Reversers', DrugSet)