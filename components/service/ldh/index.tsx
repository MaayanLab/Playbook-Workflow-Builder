import { MetaNode } from '@/spec/metanode'
import { GeneTerm } from '@/components/core/term'
import { GeneSet, RegulatoryElementSet } from '@/components/core/set'
import { RegulatoryElementSetForGeneSetInfo, MyGeneToRegulatoryElementSetInfo, MyRegulatoryElementSetInfo } from '@/components/service/regulatoryElementInfo'
import { GeneInfo, GeneInfoFromGeneTerm } from '@/components/service/mygeneinfo'
import { linkeddatahub_icon } from '@/icons'
import { z } from 'zod'
import { describe } from 'node:test'

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


function getREPositionDataFromLinkDataHub(reObjSet: any){
  let response: MyRegulatoryElementSetInfo = [];
  for(let i in reObjSet){
    let rgObj: any = reObjSet[i];

    let coordinates = null;
    if(rgObj.entContent != null && rgObj.entContent.coordinates != null){
      coordinates = rgObj.entContent.coordinates;
    }else{
      coordinates = {
          chromosome: "",
          start: 0,
          end: 0
      };
    }
    let tempObj = {
      entId: rgObj.entId,
      ldhId: rgObj.ldhId,
      entContent: {
        coordinates: coordinates
      }
    }
    response.push(tempObj);
  }
  return response;
}

export const GetRegulatoryElementsInfoForGeneInfo = MetaNode('GetRegulatoryElementsInfoForGeneInfo')
  .meta({
    label: 'Identify regulatory element in the vicinity of given gene',
    description: 'Regulatory elements in 10kbps region upstream or downstream of gene body.',
    icon: [linkeddatahub_icon],
    pagerank: 1,
  })
  .inputs({ geneInfo: GeneInfo })
  .output(RegulatoryElementSet)
  .resolve(async (props) => {
    const response =  await myGeneInfoFromLinkDataHub(props.inputs.geneInfo.symbol);
    if(response.data == null || response.data.ld == null){
      throw new Error("Unable to get data from Linked Data Hub API, please try again or wait a few minutes before the next atempt!");
    }

    let reObjSet = response.data.ld.RegulatoryElement;
    if(reObjSet == null || reObjSet.length == 0){
      throw new Error("Unable to get Regulatory Element data for gene "+props.inputs.geneInfo.symbol+" from Linked Data Hub API, please try again or wait a few minutes before the next atempt!");
    }
    
    let reIdentifiersSet = [];
    for(let i in reObjSet){
      reIdentifiersSet.push(reObjSet[i].entId);
    }
    if(reIdentifiersSet.length == 0){
      throw new Error("Unable to get Regulatory Element data for gene "+props.inputs.geneInfo.symbol+" from Linked Data Hub API, please try again or wait a few minutes before the next atempt!");
    }

    let tempObj = {
      description:"Regulatory Elements Set",
      set: reIdentifiersSet
    }
    return tempObj;
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
  }).story(props => ({
    abstract: `A list of regulatory elements in the vicinity of the gene were retrieved from the CFDE Linked Data Hub [\\ref{CFDE Linked Data Hub, https://ldh.genome.network/cfde/ldh/}].`,
    introduction: `CFDE LDH is a graph-based network that facilitates access to excerpted regulatory information from external databases and studies including SCREEN, GTEx, and EN-TEx [\\ref{doi:10.1038/s41586-020-2493-4},\\ref{doi:10.1126/science.aaz1776},\\ref{doi:10.1126/science.aar3146}].`,
    methods: `Input gene was queried through CFDE LDH API endpoints, and JSON response with linked data were parsed to retrieve FAIR identifiers for linked regulatory elements.`,
    legend: `A table displaying the globally unique identifiers and positions for regulatory elements in the gene body or in the 10kbps upstream or downstream region.`,
  })).build()

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

      let reObjSet: any = response.data.ld.RegulatoryElement;
      if(reObjSet == null || reObjSet.length == 0){
        continue;
      }
      let reInfoSet = getREPositionDataFromLinkDataHub(reObjSet);
      if(reInfoSet == null || reInfoSet.length == 0){
        continue;
      }
      let temp = {
        gene: g,
        regulatoryElements: reInfoSet
      }
      regulatoryElemForGenes.push(temp);
    }

    return regulatoryElemForGenes;
  }).story(props => ({
    abstract: `A list of regulatory elements in the vicinity of the genes were retrieved from the CFDE Linked Data Hub [\\ref{CFDE Linked Data Hub, https://ldh.genome.network/cfde/ldh/}].`,
    introduction: `CFDE LDH is a graph-based network that facilitates access to excerpted regulatory information from external databases and studies including SCREEN, GTEx, and EN-TEx [\\ref{doi:10.1038/s41586-020-2493-4},\\ref{doi:10.1126/science.aaz1776},\\ref{doi:10.1126/science.aar3146}].`,
    methods: `Input genes was queried through CFDE LDH API endpoints, and JSON response with linked data were parsed to retrieve FAIR identifiers for linked regulatory elements.`,
    legend: `A table displaying the globally unique identifiers and positions for regulatory elements in the gene body or in the 10kbps upstream or downstream region.`,
  })).build()