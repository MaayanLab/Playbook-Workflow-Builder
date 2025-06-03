import { GeneSet } from "@/components/core/set"
import { rummagene_icon } from "@/icons"
import { MetaNode } from "@/spec/metanode"
import { z } from 'zod'

const rummagene_url = 'https://rummagene.com'

export const RummageneGeneSet = MetaNode(`RummageneGeneSet`)
  .meta({
    label: 'Rummagene Gene Set',
    description: 'A gene set uploaded to rummagene',
    icon: [rummagene_icon],
    external: true,
  })
  .codec(z.object({
    id: z.string(),
  }))
  .view(props => (
    <div className="flex-grow flex flex-row m-0" style={{ minHeight: 500 }}>
      <iframe
        className="flex-grow border-0"
        src={`${rummagene_url}/enrich?dataset=${props.id}&embed`}
      />
    </div>
  ))
  .build()

export const RummageneEnrichmentAnalysis = MetaNode(`RummageneEnrichmentAnalysis`)
  .meta({
    label: 'Rummagene Enrichment Analysis',
    description: 'Use Rummagene to search through gene sets in PubMedCentral supplemental material',
    external: true,
  })
  .inputs({ gene_set: GeneSet })
  .output(RummageneGeneSet)
  .resolve(async (props) => {
    const req = await fetch(`${rummagene_url}/graphql`, {
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
      },
      method: 'POST',
      body: JSON.stringify({"operationName":"AddUserGeneSet","variables":{"description":`${props.inputs.gene_set.description ?? 'Gene set'} from pwb`,"genes":props.inputs.gene_set.set},"query":"mutation AddUserGeneSet($genes: [String], $description: String = \"\") {\n  addUserGeneSet(input: {genes: $genes, description: $description}) {\n    userGeneSet {\n      id\n      __typename\n    }\n    __typename\n  }\n}\n"}),
    })
    if (!req.ok) throw new Error('Failed to submit gene set to rummagene')
    const { data: { addUserGeneSet: { userGeneSet: { id } } } } = await req.json()
    return { id }
  })
  .story(props => ({
    abstract: `Papers in Pub Med Central with gene sets significantly overlapping with ${props.inputs?.gene_set?.description ? props.inputs.gene_set.description : 'the gene set'} were identified using Rummagene\\ref{doi:10.1038/s42003-024-06177-7}.`,
    legend: `A listing of significantly overlapping gene sets with ${props.inputs?.gene_set?.description ? props.inputs.gene_set.description : 'the gene set'} found from the supplemental material of papers using Rummagene\\ref{doi:10.1038/s42003-024-06177-7}.`,
  }))
  .build()
