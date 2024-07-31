import { MetaNode } from "@/spec/metanode";
import { EnrichrEnrichmentAnalysis } from "../enrichr";
import { view_in_graph_icon } from "@/icons";
import { z } from 'zod'

const cfde_gse_url = 'https://gse.cfde.cloud'

export const CFDEGSEKG = MetaNode('CFDEGSEKG')
  .meta({
    label: 'CFDE Gene Set Enrichment Knowledge Graph',
    description: 'A CFDE-Centric knowledge graph for gene set enrichment results',
    icon: [view_in_graph_icon],
  })
  .codec(z.union([
    z.object({
      empty: z.literal(true),
    }),
    z.object({
      shortId: z.string(),
      userListId: z.number(),
    })
  ]))
  .view(userlist => {
    const searchParams = new URLSearchParams()
    if (!('empty' in userlist)) {
      searchParams.append('q', JSON.stringify({
        "min_lib": 1,
        "libraries":[
          {"name":"LINCS_L1000_Chem_Pert_Consensus_Sigs","limit":5},
          {"name":"HuBMAP_ASCTplusB_augmented_2022","limit":5},
        ],
        "userListId":userlist.userListId.toString(),
        "search":true,
      }))
      searchParams.append('fullscreen', 'true')
    }
    return (
      <div className="flex-grow flex flex-row m-0" style={{ minHeight: 1000 }}>
        {'empty' in userlist ? 
          <div className="prose max-w-none">Enrichment Analysis Cannot be Performed on Empty Set</div>
          : <iframe
            className="flex-grow border-0"
            src={`${cfde_gse_url}?${searchParams.toString()}`}
          />}
      </div>
    )
  })
  .build()

export const EnrichrUserListToCFDEGSEKG = MetaNode('EnrichrUserListToCFDEGSEKG')
  .meta({
    label: 'Visualize CFDE Gene Set Enrichment Knowledge Graph',
    description: 'View gene set results as a knowledge graph',    
  })
  .inputs({ enrichr: EnrichrEnrichmentAnalysis })
  .output(CFDEGSEKG)
  .resolve(async (props) => props.inputs.enrichr)
  .story(props => ({
    introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
    methods: `The enrichment results from Enrichr\\ref{doi:10.1002/cpz1.90} are used to subset a Knowledge Graph of inter-connected CFDE Gene Set Libraries using CFDE GSE\\ref{CFDE GSE, https://gse.cfde.cloud/}.`,
    legend: `An interactive page containing a knowledge graph highlighting the enriched terms in the CFDE Gene Set Libraries provided by CFDE GSE\\ref{CFDE GSE, https://gse.cfde.cloud/}.`,
  }))
  .build()
