import { DrugSet,GeneSet } from "@/components/core/set"
import { ScoredGenes } from "@/components/core/scored"
import { l2s2_icon } from "@/icons"
import { MetaNode } from "@/spec/metanode"
import { z } from 'zod'

const l2s2_url = 'https://l2s2.maayanlab.cloud'

export const L2S2GeneSet = MetaNode(`L2S2GeneSet`)
  .meta({
    label: 'L2S2 Gene Set',
    description: 'A gene set uploaded to L2S2',
    icon: [l2s2_icon],
    external: true,
  })
  .codec(z.object({
    id: z.string(),
  }))
  .view(props => (
    <div className="flex-grow flex flex-row m-0" style={{ minHeight: 1000 }}>
      <iframe
        className="flex-grow border-0"
        src={`${l2s2_url}/enrich?dataset=${props.id}&embed`}
      />
    </div>
  ))
  .build()

export const L2S2GeneSignature = MetaNode(`L2S2GeneSignature`)
  .meta({
    label: 'L2S2 Gene Signature',
    description: 'A gene signature uploaded to L2S2',
    icon: [l2s2_icon],
    external: true,
  })
  .codec(z.object({ up_id: z.string(), down_id: z.string() }))
  .view(props => (
    <div className="flex-grow flex flex-row m-0" style={{ minHeight: 1000 }}>
      <iframe
        className="flex-grow border-0"
        src={`${l2s2_url}/enrichpair?dataset=${props.up_id}&dataset=${props.down_id}&embed`}
      />
    </div>
  ))
  .build()

export const L2S2EnrichmentAnalysis = MetaNode(`L2S2EnrichmentAnalysis`)
  .meta({
    label: 'LINCS L1000 Signature Search (L2S2)',
    description: 'Use L2S2 to identify small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to the gene set',
    icon: [l2s2_icon],
    external: true,
  })
  .inputs({ gene_set: GeneSet })
  .output(L2S2GeneSet)
  .resolve(async (props) => {
    const req = await fetch(`${l2s2_url}/graphql`, {
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
      },
      method: 'POST',
      body: JSON.stringify({"operationName":"AddUserGeneSet","variables":{"description":`${props.inputs.gene_set.description ?? 'Gene set'} from pwb`,"genes":props.inputs.gene_set.set},"query":"mutation AddUserGeneSet($genes: [String], $description: String = \"\") {\n  addUserGeneSet(input: {genes: $genes, description: $description}) {\n    userGeneSet {\n      id\n      __typename\n    }\n    __typename\n  }\n}\n"}),
    })
    if (!req.ok) throw new Error('Failed to submit gene set to l2s2')
    const { data: { addUserGeneSet: { userGeneSet: { id } } } } = await req.json()
    return { id }
  })
  .story(props => ({
    abstract: `Small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to gene sets with ${props.inputs?.gene_set?.description ? props.inputs.gene_set.description : 'the gene set'} were identified using L2S2\\ref{doi:10.1093/nar/gkaf373}.`,
    legend: `A listing of small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to gene sets with ${props.inputs?.gene_set?.description ? props.inputs.gene_set.description : 'the gene set'} found using L2S2\\ref{doi:10.1093/nar/gkaf373}.`,
  }))
  .build()

export const L2S2SignatureEnrichmentAnalysis = MetaNode(`L2S2SignatureEnrichmentAnalysis`)
  .meta({
    label: 'LINCS L1000 Signature Search (L2S2)',
    description: 'Use L2S2 to identify small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to the gene signature',
    icon: [l2s2_icon],
    external: true,
  })
  .inputs({ scored_genes: ScoredGenes })
  .output(L2S2GeneSignature)
  .resolve(async (props) => {
    const postSet = async (description: string, genes: string[]) => {
      const req = await fetch(`${l2s2_url}/graphql`, {
        headers: { 'Content-Type': 'application/json', 'Accept': 'application/json' },
        method: 'POST',
        body: JSON.stringify({
          operationName: 'AddUserGeneSet',
          variables: { description, genes },
          query: 'mutation AddUserGeneSet($genes: [String], $description: String = "") {\n  addUserGeneSet(input: {genes: $genes, description: $description}) {\n    userGeneSet {\n      id\n      __typename\n    }\n    __typename\n  }\n}\n',
        }),
      })
      if (!req.ok) throw new Error('Failed to submit gene set to L2S2')
      const { data: { addUserGeneSet: { userGeneSet: { id } } } } = await req.json()
      return id
    }

    const parseZ = (z:number | "nan" | "inf" | "-inf") => {
      return (z === "inf") ? Infinity : (z === "-inf") ? -Infinity : (z === "nan" || Number.isNaN(z)) ? 0 : z
    }
    const sorted = [...props.inputs.scored_genes].sort((a, b) => parseZ(b.zscore) - parseZ(a.zscore))
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
    abstract: `Small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to the signature were identified using L2S2\\ref{doi:10.1093/nar/gkaf373}.`,
    legend: `A table of enrichment results from L2S2\\ref{doi:10.1093/nar/gkaf373} showing matched LINCS L1000 signatures.`,
  }))
  .build()