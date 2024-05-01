import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/input/term'
import { z } from 'zod'
import { resolveVarinatCaID, variantIdResolveErrorMessage, alleleRegRespErrorMessage } from '../variantUtils'
import { getAlleleRegistryVariantInfo } from './alleleRegistryVariantInfo'

export const MyVariantInfoC = z.object({
    '_id':z.string()
  })
  export type MyVariantInfo = z.infer<typeof MyVariantInfoC>
  
  async function getVarinatIntoFromMyVariantInfo( link: string): Promise<MyVariantInfo> {
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

    const response = await getAlleleRegistryVariantInfo(varCaId);
    if(response == null){
      throw new Error(alleleRegRespErrorMessage);
    }

    let myVariantIdLink = null;

    if(response['externalRecords'] != null){
      let externalRecords = response['externalRecords'];
      for(let er in externalRecords){
        if(er == "MyVariantInfo_hg38"){
          myVariantIdLink = externalRecords[er][0]['@id'];
          break;
        }
      }
    }else{
      throw new Error("Unable to find requested data, missing External Records info from Allele Reg. API response!");
    }

    if(myVariantIdLink != null){
      return await getVarinatIntoFromMyVariantInfo(myVariantIdLink);
    }else{
      throw new Error("Unable to find requested data, missing MyVarintInfo API endpoint!");
    }
  })
  .story(props => `The closest gene to the variant was extract from the MyVariant.info API results [\\ref{doi:10.1093/bioinformatics/btac017}].`)
  .build()
