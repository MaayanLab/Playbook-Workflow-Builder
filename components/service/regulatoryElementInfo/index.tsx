import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { RegulatoryElementTerm } from '@/components/core/input/term'
import { GeneSet, VariantSet } from '@/components/core/input/set'
import { z } from 'zod'

export const MyRegulatoryElementC = z.object({
    data: z.object({
      entId: z.string(),
      entType: z.string(),
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

export const RegulatoryElementInfo = MetaNode.createData('RegulatoryElementInfo')
.meta({
  label: 'Regulatory Element',
  description: 'Regulatory Element resolver'
})
.codec(MyRegulatoryElementC)
.view(regElem => (
  <div> {regElem.data.entId} {regElem.data.entType} </div>  
))
.build()

export async function myRegElemInfo_query(regElemId: string): Promise<MyRegulatoryElement> {
  const res = await fetch(`https://genboree.org/cfde-gene-dev/RegulatoryElement/id/${encodeURIComponent(regElemId)}`)
  return await res.json()
}

export const RegElementInfoFromRegElementTerm = MetaNode('RegElementInfoFromRegElementTerm')
.meta({
  label: 'Resolve Regulatory Element Info from Term',
  description: 'Resolve Regultory Element info from variant term with MyVarinatInfo',
  //icon: [regulatoryElementInfo_icon],
})
.inputs({ regulatoryEelement: RegulatoryElementTerm })
.output(RegulatoryElementInfo)
.resolve(async (props) => {
  return await myRegElemInfo_query(props.inputs.regulatoryEelement);
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
  let geneInfoList = [];
  let geneList = props.inputs.regElemInfo.data.ldFor.Gene;
  for(const g of geneList){
    geneInfoList.push(g.entId);
  }
  return geneInfoList; 
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
   let variantInfoList = [];
   let geneList = props.inputs.regElemInfo.data.ld.Variant;
   for(const g of geneList){
    variantInfoList.push(g.entId);
   }
   return variantInfoList; 
})
.build()




