import { MetaNode } from '@/spec/metanode'
import { GeneTerm } from '@/components/core/input/term'
import { RegulatoryElementSet } from '@/components/core/input/set'
import { GeneInfo, GeneInfoFromGeneTerm } from '@/components/service/mygeneinfo'
import { linkeddatahub_icon } from '@/icons'
import { z } from 'zod'

export const MyGeneInfoByTermC = z.object({
  data: z.object({
    ld: z.object({
      RegulatoryElement: z.array(z.object({
            entId: z.string(),
            ldhId: z.string()
          })
        )
    })
  })
})
export type MyGeneInfoByTerm = z.infer<typeof MyGeneInfoByTermC>

export async function myGeneInfoFromLinkDataHub(geneTerm: string): Promise<MyGeneInfoByTerm> {
  const res = await fetch(`https://ldh.genome.network/cfde/ldh/Gene/id/${encodeURIComponent(geneTerm)}`)
  return await res.json()
}

export const GetRegulatoryElementsForGeneInfo = MetaNode('GetRegulatoryElementsForGeneInfo')
  .meta({
    label: 'Resolve Regulatory Elements from LDH',
    description: 'Resolve regulatory elements from gene with Linked Data Hub',
    icon: [linkeddatahub_icon],
    pagerank: 1,
  })
  .inputs({ geneInfo: GeneInfo })
  .output(RegulatoryElementSet)
  .resolve(async (props) => {
    const response =  await myGeneInfoFromLinkDataHub(props.inputs.geneInfo.symbol);
    if(response.data == null || response.data.ld == null){
      return {
        description: 'Regulatory Element set for gene is empty' ,
        set: []
      };
    }
    let reNames = response.data.ld.RegulatoryElement.map(({ entId }) => entId );
    let reSet = {
      description: 'Regulatory Element set for gene '+props.inputs.geneInfo.symbol ,
      set: reNames
    };
    return reSet;
  })
  .story(props =>
    `Regulatory elements were obtained from the Linked Data Hub [\\ref{Linked Data Hub, https://ldh.genome.network/cfde/ldh/}].`
  )
  .build()

export const GetRegulatoryElementsForGeneInfoFromGene = MetaNode('GetRegulatoryElementsForGeneInfoFromGene')
  .meta(GetRegulatoryElementsForGeneInfo.meta)
  .inputs({ gene: GeneTerm })
  .output(GetRegulatoryElementsForGeneInfo.output)
  .resolve(async (props) => {
    const geneInfo = await GeneInfoFromGeneTerm.resolve(props)
    return await GetRegulatoryElementsForGeneInfo.resolve({ ...props, inputs: { geneInfo } })
  })
  .story(props =>
    `Regulatory elements were obtained from the Linked Data Hub [\\ref{Linked Data Hub, https://ldh.genome.network/cfde/ldh/}].`
  )
  .build()
