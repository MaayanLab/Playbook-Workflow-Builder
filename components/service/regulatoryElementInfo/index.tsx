import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { RegulatoryElementTerm } from '@/components/core/input/term'
import { GeneSet, VariantSet, RegulatoryElementSet } from '@/components/core/input/set'
import { z } from 'zod'
import { linkeddatahub_icon } from '@/icons'
import { RegulatoryElementSetInfo, MyRegulatoryElementSetInfo } from '@/components/service/mygeneinfo'

export const MyRegulatoryElementC = z.object({
  data: z.object({
    entId: z.string(),
    entType: z.string(),
    coordinates: z.object({
      chromosome: z.string(),
      start: z.any(),
      end: z.any()
    }),
    ld: z.object({
      ENCODERegulatoryElementEvidence: z.array(z.object({
        ldhId: z.string()
      })),
      Variant: z.array(z.object({
        entId: z.string(),
        entIri: z.string(),
        ldhId: z.string(),
      }))
    }),
    ldFor: z.object({
      Gene: z.array(z.object({
        entId: z.string(),
        ldhId: z.string()
      }))
    })
  })
})

export type MyRegulatoryElement = z.infer<typeof MyRegulatoryElementC>

export const RE_PositionalDataC = z.object({
  data: z.object({
    cCREQuery: z.array(z.object({
      coordinates: z.object({
        chromosome: z.string(),
        start: z.number(),
        end: z.number()
      })
    }))
  })
});
type RE_PositionalData = z.infer<typeof RE_PositionalDataC>

export const RegulatoryElementInfo = MetaNode('RegulatoryElementInfo')
  .meta({
    label: 'Regulatory Element',
    description: 'Regulatory Element resolver',
    icon: [linkeddatahub_icon],
  })
  .codec(MyRegulatoryElementC)
  .view(regElem => (
    <div> {regElem.data.entId} Regulatory Element<br></br>Position: {regElem.data.coordinates.chromosome}: {regElem.data.coordinates.start}-{regElem.data.coordinates.end} (GRCh38)</div>
  ))
  .build()

export async function myRegElemInfo_query(regElemId: string): Promise<MyRegulatoryElement> {
  const res = await fetch(`https://genboree.org/cfde-gene-dev/RegulatoryElement/id/${encodeURIComponent(regElemId)}`)
  return await res.json()
}

export async function getRegElemPositionData(regElemId: string): Promise<RE_PositionalData> {
  let bodyString = '{\"query\":\"query CCRE{\\n  cCREQuery(assembly: \\"GRCh38\\", accession:\\"'+encodeURIComponent(regElemId)+'\\") {\\n  coordinates {\\n  chromosome\\n  start\\n  end\\n  }\\n  }\\n}"}';
  const res = await fetch(`https://ga.staging.wenglab.org/graphql`, {
    method: 'POST',
    headers: {
      'Accept-Encoding': 'gzip, deflate, br',
      'Content-Type': 'application/json',
      'Accept': 'application/json',
      'Connection': 'keep-alive',
      'DNT': '1',
      'Origin': 'https://ga.staging.wenglab.org'
    },
    body: bodyString
  })
  return await res.json()
}

export const GetRegElementInfoFromRegElementTerm = MetaNode('GetRegElementInfoFromRegElementTerm')
  .meta({
    label: 'Resolve Regulatory Element Info from Term',
    description: 'Resolve Regulatory Element info from variant term with MyVariantInfo',
    icon: [linkeddatahub_icon],
  })
  .inputs({ regulatoryElement: RegulatoryElementTerm })
  .output(RegulatoryElementInfo)
  .resolve(async (props) => {
    const rePositionData = await getRegElemPositionData(props.inputs.regulatoryElement);
    const response = await myRegElemInfo_query(props.inputs.regulatoryElement);
    if(rePositionData.data.cCREQuery[0].coordinates != null){
      response.data.coordinates = rePositionData.data.cCREQuery[0].coordinates;
    }else{
      response.data.coordinates = {
          chromosome: "",
          start: "",
          end: ""
      };
    }
    return response;
  })
  .build()

export const GetRegElementSetInfoFromRegElementTerm = MetaNode('GetRegElementSetInfoFromRegElementTerm')
  .meta({
    label: 'Resolve Regulatory Element Set Info from Term',
    description: 'Resolve Regulatory Element Set Info from term (id).',
    icon: [linkeddatahub_icon],
  })
  .inputs({ regulatoryElementSet: RegulatoryElementSet })
  .output(RegulatoryElementSetInfo)
  .resolve(async (props) => {
    let regElemeIdsSet = props.inputs.regulatoryElementSet.set;
    let response: MyRegulatoryElementSetInfo = [];
    for(let i in regElemeIdsSet){
      let rgId = regElemeIdsSet[i];
      const rePositionData = await getRegElemPositionData(rgId);
      let coordinates = null;
      if(rePositionData.data.cCREQuery[0] != null && rePositionData.data.cCREQuery[0].coordinates != null){
        coordinates = rePositionData.data.cCREQuery[0].coordinates;
      }else{
        coordinates = {
            chromosome: "",
            start: 0,
            end: 0
        };
      }
      let rgObj = {
        entId: rgId,
        ldhId: rgId,
        entContent: {
          coordinates: coordinates
        }
      }
      response.push(rgObj);
    }

    return response;
  })
  .build()

export const GetGenesForRegulatoryElementInfo = MetaNode('GetGenesForRegulatoryElementInfo')
  .meta({
    label: 'Resolve Genes for Regulatory Elements Info',
    description: 'Get Linked Genes For Regulatory Element Info',
  })
  .inputs({ regElemInfo: RegulatoryElementInfo  })
  .output(GeneSet)
  .resolve(async (props) => {
    let geneNames =  props.inputs.regElemInfo.data.ldFor.Gene.map(({ entId }) => entId);
    let geneSet = {
      description: 'Gene set for regulatory element '+props.inputs.regElemInfo.data.entId,
      set: geneNames
    };
    return geneSet;
  })
  .build()

export const GetVariantsForRegulatoryElementInfo = MetaNode('GetVariantListForRegulatoryElementInfo')
  .meta({
    label: 'Resolve Variants for Regulatory Elements Info',
    description: 'Get Linked Variants For Regulatory Element Info',
  })
  .inputs({ regElemInfo: RegulatoryElementInfo  })
  .output(VariantSet)
  .resolve(async (props) => {
    let variantNames =  props.inputs.regElemInfo.data.ld.Variant.map(({ entId }) => entId);
    let variantSet = {
      description: 'Variant set for regulatory element '+props.inputs.regElemInfo.data.entId,
      set: variantNames
    };
    return variantSet;
  })
  .build()
