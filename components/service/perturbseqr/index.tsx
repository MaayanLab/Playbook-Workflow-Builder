import { DrugSet, GeneSet } from "@/components/core/set"
import { ScoredGenes } from "@/components/core/scored"
import python from '@/utils/python'
import { perturbseqr_icon } from "@/icons"
import { MetaNode } from "@/spec/metanode"
import { z } from 'zod'

const perturbseqr_url = 'https://perturbseqr.maayanlab.cloud'

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
    <div className="flex-grow flex flex-row m-0" style={{ minHeight: 1000 }}>
      <iframe
        className="flex-grow border-0"
        src={`${perturbseqr_url}/enrich?dataset=${props.id}&embed`}
      />
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
    <div className="flex-grow flex flex-row m-0" style={{ minHeight: 1000 }}>
      <iframe
        className="flex-grow border-0"
        src={`${perturbseqr_url}/enrichpair?dataset=${props.up_id}&dataset=${props.down_id}&embed`}
      />
    </div>
  ))
  .build()

export const PerturbSeqrEnrichmentAnalysis = MetaNode(`Perturb-SeqrEnrichmentAnalysis`)
  .meta({
    label: 'Perturb-Seqr Gene Set Search',
    description: 'Use Perturb-Seqr to identify small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to the gene set',
    icon: [perturbseqr_icon],
    external: true,
  })
  .inputs({ gene_set: GeneSet })
  .output(PerturbSeqrGeneSet)
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
    return { id }
  })
  .story(props => ({
    abstract: `Small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to gene sets with ${props.inputs?.gene_set?.description ? props.inputs.gene_set.description : 'the gene set'} were identified using Perturb-Seqr\\ref{Perturb-Seqr, https://perturbseqr.maayanlab.cloud/}.`,
    legend: `A listing of small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to gene sets with ${props.inputs?.gene_set?.description ? props.inputs.gene_set.description : 'the gene set'} found using Perturb-Seqr\\ref{Perturb-Seqr, https://perturbseqr.maayanlab.cloud/}.`,
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
  .output(PerturbSeqrGeneSignature)
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

    const sorted = [...props.inputs.scored_genes].sort((a, b) => +b.zscore - +a.zscore)
    const [up_id, down_id] = await Promise.all([
      postSet('Up gene set from pwb', sorted
        .filter(g => g.zscore === 'inf' || (typeof g.zscore === 'number' && g.zscore > 0))
        .map(g => g.term)
        .slice(0,500)),
      postSet('Down gene set from pwb', sorted
        .filter(g => g.zscore === '-inf' || (typeof g.zscore === 'number' && g.zscore < 0))
        .map(g => g.term)
        .slice(-500)),
    ])

    return { up_id:up_id, down_id:down_id }
  })
  .story(props => ({
    abstract: `Small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to the signature were identified using Perturb-Seqr\\ref{Perturb-Seqr, https://perturbseqr.maayanlab.cloud/}.`,
    legend: `A table of drug and gene perturbation enrichment results from Perturb-Seqr\\ref{Perturb-Seqr, https://perturbseqr.maayanlab.cloud/}.`,
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
    .inputs({ signature: PerturbSeqrGeneSignature })
    .output(SetT)
    .resolve(async (props) => {
      return await python(
        'components.service.perturbseqr.extract_perturbation_set',
        { kargs: [props.inputs.signature, perturbationType, direction] },
        message => props.notify({ type: 'info', message }),
      )
    })
    .story(props => ({
      abstract: `${perturbationType} ${direction.toLocaleLowerCase()} were extracted from the Perturb-Seqr\\ref{Perturb-Seqr, https://perturbseqr.maayanlab.cloud/} enrichment results.`,
      methods: `Rows whose dataset matched the ${perturbationType.toLowerCase()} perturbation libraries with p-value < 0.05 were selected from the Perturb-Seqr\\ref{Perturb-Seqr, https://perturbseqr.maayanlab.cloud/} enrichment results.`,
      tableLegend: `A set of ${perturbationType.toLocaleLowerCase()} ${direction.toLocaleLowerCase()} identified from the Perturb-Seqr\\ref{Perturb-Seqr, https://perturbseqr.maayanlab.cloud/} enrichment results.`,
    }))
    .build()
}

export const ExtractPerturbSeqrMimickerGenes = makePerturbSeqrExtractNode('Gene', 'Mimickers', GeneSet)
export const ExtractPerturbSeqrMimickerDrugs = makePerturbSeqrExtractNode('Drug', 'Mimickers', DrugSet)
export const ExtractPerturbSeqrReverserGenes = makePerturbSeqrExtractNode('Gene', 'Reversers', GeneSet)
export const ExtractPerturbSeqrReverserDrugs = makePerturbSeqrExtractNode('Drug', 'Reversers', DrugSet)