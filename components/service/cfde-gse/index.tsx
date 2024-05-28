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
  .story(props => ``)
  .build()
