import { GeneSet } from "@/components/core/set"
import { ScoredGenes } from "@/components/core/scored"
import dynamic from 'next/dynamic'
import python from '@/utils/python'
import { downloadBlob } from '@/utils/download'
import { l2s2_icon } from "@/icons"
import { MetaNode } from "@/spec/metanode"
import { z } from 'zod'

const l2s2_url = 'https://l2s2.maayanlab.cloud'
const Matrix = dynamic(() => import('@/app/components/Matrix'))

const L2S2Value = z.union([z.string(), z.number(), z.boolean(), z.null()])

export const L2S2EnrichmentResults = MetaNode(`L2S2EnrichmentResults`)
  .meta({
    label: 'L2S2 Enrichment Results',
    description: 'Enrichment results from L2S2',
    icon: [l2s2_icon],
    external: true,
  })
  .codec(z.object({
    shape: z.tuple([z.number(), z.number()]),
    index: z.array(z.string()),
    columns: z.array(z.string()),
    values: z.array(z.array(L2S2Value)),
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
          'l2s2_enrichment.json'
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
          'l2s2_enrichment.csv'
        ),
      }}
    />
  ))
  
  .build()

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
    <div className="prose max-w-none">
      <p>Gene set uploaded to L2S2 with ID: <code>{props.id}</code></p>
      <p>Use <em>Fetch L2S2 Enrichment Results</em> to load the results table.</p>
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
    <div className="prose max-w-none">
      <p>Gene signature uploaded to L2S2.</p>
      <p>Up ID: <code>{props.up_id}</code> / Down ID: <code>{props.down_id}</code></p>
      <p>Use <em>Fetch L2S2 Enrichment Results</em> to load the results table.</p>
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
  .output(L2S2EnrichmentResults)
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
    return await python(
      'components.service.l2s2.fetch_geneset_enrichment_results',
      { kargs: [id] },
      message => props.notify({ type: 'info', message }),
    )
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
  .output(L2S2EnrichmentResults)
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

    const [up_id, down_id] = await Promise.all([
      postSet('Up gene set from pwb', props.inputs.scored_genes
        .filter(g => g.zscore === 'inf' || (typeof g.zscore === 'number' && g.zscore > 0))
        .map(g => g.term)),
      postSet('Down gene set from pwb', props.inputs.scored_genes
        .filter(g => g.zscore === '-inf' || (typeof g.zscore === 'number' && g.zscore < 0))
        .map(g => g.term)),
    ])

    return await python(
      'components.service.l2s2.fetch_enrichment_results',
      { kargs: [up_id, down_id] },
      message => props.notify({ type: 'info', message }),
    )
  })
  .story(props => ({
    abstract: `Small molecules and single gene CRISPR KOs that produce gene expression profiles similar or opposite to the signature were identified using L2S2\\ref{doi:10.1093/nar/gkaf373}.`,
    legend: `A table of enrichment results from L2S2\\ref{doi:10.1093/nar/gkaf373} showing matched LINCS L1000 signatures.`,
  }))
  .build()