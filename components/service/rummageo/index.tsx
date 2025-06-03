import { GeneSet } from "@/components/core/set"
import { rummageo_icon } from "@/icons"
import { MetaNode } from "@/spec/metanode"
import { z } from 'zod'

const rummageo_url = 'https://rummageo.com'

export const RummaGEOGeneSet = MetaNode(`RummaGEOGeneSet`)
  .meta({
    label: 'RummaGEO Gene Set',
    description: 'A gene set uploaded to RummaGEO',
    icon: [rummageo_icon],
    external: true,
  })
  .codec(z.object({
    id: z.string(),
  }))
  .view(props => (
    <div className="flex-grow flex flex-row m-0" style={{ minHeight: 1000 }}>
      <iframe
        className="flex-grow border-0"
        src={`${rummageo_url}/enrich?dataset=${props.id}&embed`}
      />
    </div>
  ))
  .build()

export const RummaGEOEnrichmentAnalysis = MetaNode(`RummaGEOEnrichmentAnalysis`)
  .meta({
    label: 'RummaGEO Enrichment Analysis',
    description: 'Use RummaGEO to search through differentially expressed GEO study gene sets',
    external: true,
  })
  .inputs({ gene_set: GeneSet })
  .output(RummaGEOGeneSet)
  .resolve(async (props) => {
    const req = await fetch(`${rummageo_url}/graphql`, {
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
      },
      method: 'POST',
      body: JSON.stringify({"operationName":"AddUserGeneSet","variables":{"description":`${props.inputs.gene_set.description ?? 'Gene set'} from pwb`,"genes":props.inputs.gene_set.set},"query":"mutation AddUserGeneSet($genes: [String], $description: String = \"\") {\n  addUserGeneSet(input: {genes: $genes, description: $description}) {\n    userGeneSet {\n      id\n      __typename\n    }\n    __typename\n  }\n}\n"}),
    })
    if (!req.ok) throw new Error('Failed to submit gene set to RummaGEO')
    const { data: { addUserGeneSet: { userGeneSet: { id } } } } = await req.json()
    return { id }
  })
  .story(props => ({
    abstract: `Automatically computed signatures from GEO studies significantly overlapping with ${props.inputs?.gene_set?.description ? props.inputs.gene_set.description : 'the gene set'} were identified using RummaGEO\\ref{doi:10.1016/j.patter.2024.101072}.`,
    legend: `A listing of significantly overlapping gene sets with ${props.inputs?.gene_set?.description ? props.inputs.gene_set.description : 'the gene set'} found from differentially expressed genes from GEO studies using RummaGEO\\ref{doi:10.1016/j.patter.2024.101072}.`,
  }))
  .build()
