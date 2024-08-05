import { z } from 'zod'
import { alleleRegRespErrorMessage } from '../variantUtils'

const AlleleRegistryVariantInfoC = z.object({
    '@id': z.string(),
    communityStandardTitle: z.array(z.string()),
    entId: z.string(),
    externalRecords: z.any({}).optional()
  })
export type AlleleRegistryVariantInfo = z.infer<typeof AlleleRegistryVariantInfoC>

export async function getAlleleRegistryVariantInfo(variantId: string): Promise<AlleleRegistryVariantInfo> {
  const res = await fetch(`https://reg.genome.network/allele/${encodeURIComponent(variantId)}`);
  return await res.json();
}

export async function getVariantSetInfo(variantSetInpt: String[]){
  if(variantSetInpt == null || variantSetInpt.length == 0){
    return;
  }
  
  var varInfoSet = [];
  var varAlleleRegResponse = null;
  for(let indx in variantSetInpt){
    varAlleleRegResponse = await getAlleleRegistryVariantInfo(variantSetInpt[indx].trim());
    if(varAlleleRegResponse == null){
      continue;
    }
    varAlleleRegResponse.entId = variantSetInpt[indx].trim();
    varInfoSet.push(varAlleleRegResponse);
  }

  if(varInfoSet.length == 0){
    throw new Error(alleleRegRespErrorMessage);
  }
  return varInfoSet;
}
