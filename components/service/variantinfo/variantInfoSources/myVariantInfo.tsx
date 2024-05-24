import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/term'
import { z } from 'zod'
import { resolveVarinatCaID, variantIdResolveErrorMessage, alleleRegRespErrorMessage, getMyVarintInfoLink } from '../variantUtils'
import { getAlleleRegistryVariantInfo } from './alleleRegistryVariantInfo'

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
  export const VariantInfoFromVariantTermMyVarintInfo = MetaNode('VariantInfoFromVariantTermMyVarintInfo')
  .meta({
    label: 'Resolve Variant Info from Term (MyVarintInfo)',
    description: 'Resolve variant info from variant term using the MyVarintInfo API.',
  })
  .inputs({ variant: VariantTerm })
  .output(MyVariantInfo)
  .resolve(async (props) => {
    var varCaId = await resolveVarinatCaID(props.inputs.variant);
    if(varCaId == null || varCaId == ''){
      throw new Error(variantIdResolveErrorMessage);
    }

    const alleleRegResponse = await getAlleleRegistryVariantInfo(varCaId);
    if(alleleRegResponse == null){
      throw new Error(alleleRegRespErrorMessage);
    }

    let myVariantInfoURL = getMyVarintInfoLink(alleleRegResponse);

    if(myVariantInfoURL != null){
      return await getVarinatInfoFromMyVariantInfo(myVariantInfoURL);
    }else{
      throw new Error("Unable to find requested data, missing MyVarintInfo API link in External resources!");
    }
  })
  .story(props => ``)
  .build()
