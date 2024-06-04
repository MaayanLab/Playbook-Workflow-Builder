import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'

  export const MyVariantInfoC = z.object({
    '_id':z.string()
  })
  export type MyVariantInfo = z.infer<typeof MyVariantInfoC>
  
  export function assembleMyVarintInfoLinkWithHGVS(hgvs: string){
    return "http://myvariant.info/v1/variant/"+hgvs+"?assembly=hg38";
  }

  export async function getVariantInfoFromMyVariantInfo( link: string): Promise<MyVariantInfo> {
    const req = await fetch(link);
    return await req.json()
  }
  export const MyVariantInfo = MetaNode('MyVariantInfo')
  .meta({
    label: 'Variant Information',
    description: 'A Variant resolved with MyVariantInfo'
  })
  .codec(MyVariantInfoC)
  .view(variantinfo => (
    <div className="prose max-w-none">
      {variantinfo._id}
    </div>
  ))
  .build()
