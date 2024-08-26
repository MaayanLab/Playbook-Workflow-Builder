import { MetaNode } from '@/spec/metanode'
import { GeneTerm } from '@/components/core/term'
import { GeneSet } from '@/components/core/set'
import { RegulatoryElementSetInfo, RegulatoryElementSetForGeneSetInfo, setGenomicPositionsForRegulatoryElementSet, MyGeneToRegulatoryElementSetInfo } from '@/components/service/regulatoryElementInfo'
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

export const GetRegulatoryElementsInfoForGeneInfo = MetaNode('GetRegulatoryElementsInfoForGeneInfo')
  .meta({
    label: 'Identify regulatory element in the vicinity of given gene',
    description: 'Regulatory elements in 10kbps region upstream or downstream of gene body.',
    icon: [linkeddatahub_icon],
    pagerank: 1,
  })
  .inputs({ geneInfo: GeneInfo })
  .output(RegulatoryElementSetInfo)
  .resolve(async (props) => {
    const response =  await myGeneInfoFromLinkDataHub(props.inputs.geneInfo.symbol);
    if(response.data == null || response.data.ld == null){
      throw new Error("Unable to get data from Linked Data Hub API, please try again or wait a few minutes before the next atempt!");
    }
    let reNamesSet: string[] = response.data.ld.RegulatoryElement.map(({ entId }) => entId );
    return await setGenomicPositionsForRegulatoryElementSet(reNamesSet);
  })
  .story(props => ({
    abstract: `Regulatory elements in 10kbps region upstream or downstream of gene ${props.inputs ? ` ${props.inputs.geneInfo.symbol}` : ''}.`
  }))
  .build()

  export const GetRegulatoryElementsForGeneInfoFromGene = MetaNode('GetRegulatoryElementsForGeneInfoFromGene')
  .meta(GetRegulatoryElementsInfoForGeneInfo.meta)
  .inputs({ gene: GeneTerm })
  .output(GetRegulatoryElementsInfoForGeneInfo.output)
  .resolve(async (props) => {
    const geneInfo = await GeneInfoFromGeneTerm.resolve(props)
    return await GetRegulatoryElementsInfoForGeneInfo.resolve({ ...props, inputs: { geneInfo } })
  })
  .story(props => ({
    abstract: `Regulatory elements were obtained from the CFDE Linked Data Hub\\ref{CFDE Linked Data Hub, https://ldh.genome.network/cfde/ldh/}.`
  }))
  .build()

  export const GetRegulatoryElementsInfoForGeneInfoSet = MetaNode('GetRegulatoryElementsInfoForGeneInfoSet')
  .meta({
    label: 'Identify regulatory elements in the vicinity of given gene(s)',
    description: 'Regulatory elements in 10kbps region upstream or downstream of gene body.',
    icon: [linkeddatahub_icon],
    pagerank: 1,
  })
  .inputs({ geneSet: GeneSet })
  .output(RegulatoryElementSetForGeneSetInfo)
  .resolve(async (props) => {
    let geteSet = props.inputs.geneSet.set;

    let regulatoryElemForGenes: MyGeneToRegulatoryElementSetInfo = []; 
    for(let indx in geteSet){
      let g = geteSet[indx];

      const response =  await myGeneInfoFromLinkDataHub(g);
      if(response.data == null || response.data.ld == null){
        continue;
      }

      let reNamesSet: string[] = response.data.ld.RegulatoryElement.map(({ entId }) => entId );
      let reInfoSet = await setGenomicPositionsForRegulatoryElementSet(reNamesSet);

      let temp = {
        gene: g,
        regulatoryElements: reInfoSet
      }
      regulatoryElemForGenes.push(temp);
    }

    return regulatoryElemForGenes;
  })
  .story(props => ({
    abstract: `Regulatory elements in 10kbps region upstream or downstream of genes.`
  }))
  .build()
