export let variantIdResolveErrorMessage = "Unable to resolve the suplied Variant ID to a CAID!";
export let alleleRegRespErrorMessage = "Unable to get data from the Allele Registry API, please try again or wait a few minutes before the next attempt!";
export let gitDataHubErroMessage = "Unable to get data from Git Data Hub API, please try again or wait a few minutes before the next attempt!";


let caIdRegex = "^(CA|ca)[0-9]";
let rsIdRegex = "^(RS|rs)[0-9]";
let hgvsRegex = "^NM_|XM_[0-9]\.[0-9]:$";
let clinvarRegex = "^[0-9]+$";
let gnomADRegex = "[0-9]*-[0-9]*-[A-Za-z]*-[A-Za-z]*";
let myVariantInfoHG38Regex = "chr[0-9]*:[a-z]\.[0-9]*[A-Za-z]*";

function getCaIdFromAlleleRegistryLink(jsonObj: any){
    let alleleRegistrylink = null;
    if(Array.isArray(jsonObj)){
      alleleRegistrylink = (jsonObj[0])['@id'];
    }else{
      alleleRegistrylink = jsonObj['@id'];
    }
    let alleleRegistrylinkArray = alleleRegistrylink.split("/")
    return alleleRegistrylinkArray[(alleleRegistrylinkArray.length)-1];
  }
  
  async function getGenomeNetworkAlleles_dbSNP(variantIdTerm: string){
    const res = await fetch(`https://reg.test.genome.network/alleles?dbSNP.rs=${encodeURIComponent(variantIdTerm)}`)
    let respObj = await res.json();
    return getCaIdFromAlleleRegistryLink(respObj);
  }
  
  async function getGenomeNetworkAlleles_HGVS(variantIdTerm: string){
    const res = await fetch(`https://reg.test.genome.network/allele?hgvs=${encodeURIComponent(variantIdTerm)}`)
    let respObj = await res.json();
    return getCaIdFromAlleleRegistryLink(respObj);
  }
  
  async function getGenomeNetworkAlleles_ClinVar(variantIdTerm: string){
    const res = await fetch(`https://reg.test.genome.network/alleles?ClinVar.variationId=${encodeURIComponent(variantIdTerm)}`)
    let respObj = await res.json();
    return getCaIdFromAlleleRegistryLink(respObj);
  }
  
  async function getGenomeNetworkAlleles_gnomAD(variantIdTerm: string){
    const res = await fetch(`https://reg.test.genome.network/alleles?gnomAD.id=${encodeURIComponent(variantIdTerm)}`)
    let respObj = await res.json();
    return getCaIdFromAlleleRegistryLink(respObj);
  }
  
  async function getGenomeNetworkAlleles_myVariantInfoHG38(variantIdTerm: string){
    const res = await fetch(`https://reg.test.genome.network/alleles?MyVariantInfo_hg38.id=${encodeURIComponent(variantIdTerm)}`)
    let respObj = await res.json();
    return getCaIdFromAlleleRegistryLink(respObj);
  }
  
  export async function resolveVarinatCaID(variantIdTerm: string){
    let caIdRegexObj = new RegExp(caIdRegex);
    let rsIdRegexObj = new RegExp(rsIdRegex);
    let hgvsRegexObj = new RegExp(hgvsRegex);
    let clinvarRegexObj = new RegExp(clinvarRegex);
    let gnomADRegexObj = new RegExp(gnomADRegex);
    let myVariantInfoHG38RegexObj = new RegExp(myVariantInfoHG38Regex);
  
    if(caIdRegexObj.test(variantIdTerm)){
      return variantIdTerm;
    }else if(rsIdRegexObj.test(variantIdTerm)){
      return await getGenomeNetworkAlleles_dbSNP(variantIdTerm.slice(2));
    }else if(hgvsRegexObj.test(variantIdTerm)){
      return await getGenomeNetworkAlleles_HGVS(variantIdTerm);
    }else if(clinvarRegexObj.test(variantIdTerm)){
      return await getGenomeNetworkAlleles_ClinVar(variantIdTerm);
    }else if(gnomADRegexObj.test(variantIdTerm)){
      return await getGenomeNetworkAlleles_gnomAD(variantIdTerm);
    }else if(myVariantInfoHG38RegexObj.test(variantIdTerm)){
      return await getGenomeNetworkAlleles_myVariantInfoHG38(variantIdTerm);
    }else{
      return variantIdTerm;
    }
  }

  export function validateHGVSInput(variantIdTerm: string){
    let hgvsRegexObj = new RegExp(hgvsRegex);

    if(hgvsRegexObj.test(variantIdTerm)){
      return true;
    }else{
      return false;
    }
  }

  export function getMyVarintInfoLink(alleleRegResponse: any){
    let myVariantIdLink = null;
    if(alleleRegResponse['externalRecords'] != null){
      let externalRecords = alleleRegResponse['externalRecords'];
      for(let er in externalRecords){
        if(er == "MyVariantInfo_hg38"){
          myVariantIdLink = externalRecords[er][0]['@id'];
          break;
        }
      }
    }else{
      throw new Error("Unable to find requested data, missing External Records info from Allele Reg. API response!");
    }

    return myVariantIdLink;
  }
