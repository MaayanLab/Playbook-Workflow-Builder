import { z } from 'zod'
import { alleleRegRespErrorMessage } from '../variantUtils'

const AlleleRegistryVariantInfoC = z.object({
    '@id': z.string(),
    entId: z.string(),
    externalRecords: z.any({}).optional()
  })
export type AlleleRegistryVariantInfo = z.infer<typeof AlleleRegistryVariantInfoC>

export async function getAlleleRegistryVariantInfo(variantId: string): Promise<AlleleRegistryVariantInfo> {
  const res = await fetch(`https://reg.genome.network/allele/${encodeURIComponent(variantId)}`);
  return await res.json();
}

export async function getVariantSetInfo(varinatSetInpt: String[]){
  if(varinatSetInpt == null || varinatSetInpt.length == 0){
    return;
  }
  
  var varInfoSet = [];
  var varAlleleRegResponse = null;
  for(let indx in varinatSetInpt){
    varAlleleRegResponse = await getAlleleRegistryVariantInfo(varinatSetInpt[indx].trim());
    if(varAlleleRegResponse == null){
      continue;
    }
    varAlleleRegResponse.entId = varinatSetInpt[indx].trim();
    varInfoSet.push(varAlleleRegResponse);
  }

  if(varInfoSet.length == 0){
    throw new Error(alleleRegRespErrorMessage);
  }
  return varInfoSet;
}
