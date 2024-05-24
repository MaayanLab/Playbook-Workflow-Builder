import { z } from 'zod'

  export const MyVariantInfoC = z.object({
    '_id':z.string()
  })
  export type MyVariantInfo = z.infer<typeof MyVariantInfoC>
  
  export function assembleMyVarintInfoLinkWithHGVS(hgvs: string){
    return "http://myvariant.info/v1/variant/"+hgvs+"?assembly=hg38";
  }

  export async function getVarinatInfoFromMyVariantInfo( link: string): Promise<MyVariantInfo> {
    const req = await fetch(link);
    return await req.json()
  }
