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
            z.object({ entId: z.string() })
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
