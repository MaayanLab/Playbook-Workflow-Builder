import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/input/term'
import { xQTL_EvidenceSet } from '@/components/core/input/set'
import { z } from 'zod'
//import { varinat_icon, variantinfo_icon } from '@/icons'


export const MyVariantInfoC =  z.object({
  data: z.object({
    entId: z.string(),
    entType: z.string(),
    ld: z.object({
      xqtlEvidence: z.array(
        z.object({   
          entId: z.string(), 
          ldhId: z.string(), 
          ldhIri: z.string(),
          entContent: z.object({
            GTExIri: z.string(),
            score: z.number(),
            sourceDescription: z.string()
            /*
            eQTL: z.object({
              nes: z.string(),
              sig: z.string()
            }),
            sQTL: z.object({
              nes: z.string(),
              sig: z.string()
            })*/
          })
        })
      )
    }),
    ldFor: z.object({
      RegulatoryElement: z.array(
          z.object({ entId: z.string() })
        )
    })
  })
})

export type MyVariantInfo = z.infer<typeof MyVariantInfoC>

export const VariantInfo = MetaNode('VariantInfo')
  .meta({
    label: 'Variant Iformation',
    description: 'A Variant resolved with MyVariantInfo',
    //icon: [varinat_icon]
  })
  .codec(MyVariantInfoC)
  .view(varinatinfo => (
    <div>
      <a target="_blank" href={`https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=${varinatinfo.data.entId}`}>{varinatinfo.data.entId}</a> (variant)
    </div>  
  ))
  .build()

  async function myvariantinfo_query(variantId: string): Promise<MyVariantInfo> {
    const res = await fetch(`https://genboree.org/cfde-gene-dev/Variant/id/${encodeURIComponent(variantId)}`)
    return await res.json()
  }

  export const VarinatInfoFromVariantTerm = MetaNode('VarinatInfoFromVariantTerm')
  .meta({
    label: 'Resolve Variant Info from Term',
    description: 'Resolve variant info from variant term with MyVarinatInfo',
    //icon: [variantinfo_icon],
  })
  .inputs({ variant: VariantTerm })
  .output(VariantInfo)
  .resolve(async (props) => {
    return await myvariantinfo_query(props.inputs.variant);
  })
  .build()

  export const GetxQTL_EvidencesSetForVariantInfo = MetaNode('GetxQTL_EvidencesSetForVariantInfo')
  .meta({
    label: 'Resolve xQTL Evidence Set for Variant Info',
    description: 'Resolve xQTL Evidence for Variant Info Data',
  })
  .inputs({ variantInfo: VariantInfo  })
  .output(xQTL_EvidenceSet)
  .resolve(async (props) => {
    return props.inputs.variantInfo.data.ld.xqtlEvidence.map(({ entId }) => entId);
  })
  .build()

  export const xQTL_EvidenceDataTable = MetaNode('xQTL_EvidenceDataTable')
  .meta({
    label: 'xQTL_EvidenceDataTable',
    description: ''
  })
  .codec(MyVariantInfoC)
  .view(varinatinfo => (
    <table style={{borderCollapse: 'collapse'}}>
    <thead>
      <tr>
        <th>LDH Id</th>
        <th>GTEx Portal</th>
        <th>Score</th>
        <th>Source Description</th>     
      </tr>
    </thead>
    <tbody>
      {varinatinfo.data.ld.xqtlEvidence.map(xqtlEvidence =>
        <tr style={{borderTop: '1px solid lightgrey'}}>
          <td style={{verticalAlign: 'center', borderRight: '2px solid lightgrey', borderLeft: '2px solid lightgrey', minWidth: '100px'}}>{xqtlEvidence.ldhId}</td>
          <td style={{verticalAlign: 'center', borderRight: '2px solid lightgrey', minWidth: '100px'}}><a target="_blank" href={`${xqtlEvidence.entContent.GTExIri}`}>evidence link</a></td>
          <td style={{verticalAlign: 'center', borderRight: '2px solid lightgrey', minWidth: '70px', fontWeight:'normal'}}>{xqtlEvidence.entContent.score}</td>
          <td style={{verticalAlign: 'center', borderRight: '2px solid lightgrey', minWidth: '100px', fontWeight:'normal'}}>{xqtlEvidence.entContent.sourceDescription}</td>
        </tr>
      )}
    </tbody>
    </table>
  ))
  .build()
  /*
  <th>Normalized Effect Size (nes)</th>
  <th>p-value (sig)</th>
  <td style={{verticalAlign: 'center', borderRight: '2px solid lightgrey', minWidth: '70px', fontWeight:'normal'}}>{xqtlEvidence.entContent.eQTL.nes}</td>
  <td style={{verticalAlign: 'center', borderRight: '2px solid lightgrey', minWidth: '70px', fontWeight:'normal'}}>{xqtlEvidence.entContent.eQTL.sig}</td>
  */

  export const GetxQTL_EvidencesDataForVariantInfo = MetaNode('GetxQTL_EvidencesDataForVariantInfo')
  .meta({
    label: 'Resolve xQTL Evidence Data for Variant Info',
    description: 'Resolve xQTL Evidence Data for Variant Info Data',
  })
  .inputs({ variantInfo: VariantInfo  })
  .output(xQTL_EvidenceDataTable)
  .resolve(async (props) => {
    return props.inputs.variantInfo;
  })
  .build()
