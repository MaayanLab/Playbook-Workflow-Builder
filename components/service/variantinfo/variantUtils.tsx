export let variantIdResolveErrorMessage = "Unable to resolve the inputted Variant ID to a CAID!";
export let alleleRegRespErrorMessage = "Unable to get data from the Allele Registry API, please try again or wait a few minutes before the next attempt!";
export let linkedDataHubErroMessage = "Unable to get data from Linked Data Hub API or the variant might not be in CFDE LDH, please check the inputted identifier and try again or wait a few minutes before the next attempt!";

let variantIdentifierBaseURL = "https://reg.test.genome.network/";

//example: CA321211, CA12345
let caIdRegex = "^(CA|ca)[0-9]"; 
let caIdRegexObj = new RegExp(caIdRegex);

//example: RS369602258 
let rsIdRegex = "^(RS|rs)[0-9]";
let rsIdRegexObj = new RegExp(rsIdRegex);

//example: NM_002496.3:c.64C>T
let hgvsRegex = "^NM_|XM_[0-9]\.[0-9]:$";
let hgvsRegexObj = new RegExp(hgvsRegex);

//example: 214835
let clinvarRegex = "^[0-9]+$";
let clinvarRegexObj = new RegExp(clinvarRegex);

//example: RCV000276295
let clinvarRCVRegex = "^(RCV|rcv)[0-9]"; 
let clinvarRCVRegexObj = new RegExp(clinvarRCVRegex);

//example: CACN350004729
let cacnRegex = "^(CACN|cacn)[0-9]"; 
let cacnRegexObj = new RegExp(cacnRegex);

//example: 5-112043382-A-T
let gnomADRegex = "[0-9]*-[0-9]*-[A-Za-z]*-[A-Za-z]*";
let gnomADRegexObj = new RegExp(gnomADRegex);

//example: chr9:g.107620835G>A
let myVariantInfoHG38Regex = "chr[0-9]*:[a-z]\.[0-9]*[A-Za-z]*";
let myVariantInfoHG38RegexObj = new RegExp(myVariantInfoHG38Regex);

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

  async function getGenomeNetworkAllelesForSpecificIdentifier(variantIdTerm: string, urlProperty: string){
    const res = await fetch(variantIdentifierBaseURL+urlProperty+`=${encodeURIComponent(variantIdTerm)}`)
    let respObj = await res.json();
    return getCaIdFromAlleleRegistryLink(respObj);
  }

  export async function resolveVariantCaID(variantIdTerm: string){
    if(caIdRegexObj.test(variantIdTerm)){
      return variantIdTerm;
    }else if(rsIdRegexObj.test(variantIdTerm)){
      return await getGenomeNetworkAllelesForSpecificIdentifier(variantIdTerm.slice(2), 'alleles?dbSNP.rs');
    }else if(hgvsRegexObj.test(variantIdTerm)){
      return await getGenomeNetworkAllelesForSpecificIdentifier(variantIdTerm, 'allele?hgvs');
    }else if(clinvarRegexObj.test(variantIdTerm)){
      return await  getGenomeNetworkAllelesForSpecificIdentifier(variantIdTerm, 'alleles?ClinVar.variationId');
    }else if(clinvarRCVRegexObj.test(variantIdTerm)){
      return await getGenomeNetworkAllelesForSpecificIdentifier(variantIdTerm, 'alleles?ClinVar.RCV');
    }else if(cacnRegexObj.test(variantIdTerm)){
      return await getGenomeNetworkAllelesForSpecificIdentifier(variantIdTerm, 'by_cacnid?cacnid');
    }else if(gnomADRegexObj.test(variantIdTerm)){
      return await getGenomeNetworkAllelesForSpecificIdentifier(variantIdTerm, 'alleles?gnomAD.id');
    }else if(myVariantInfoHG38RegexObj.test(variantIdTerm)){
      return await getGenomeNetworkAllelesForSpecificIdentifier(variantIdTerm, 'alleles?MyVariantInfo_hg38.id');
    }else{
      return null;
    }
  }

  export function validateMyVariantInfoInput(variantIdTerm: string){
    let myVariantInfoHG38RegexObj = new RegExp(myVariantInfoHG38Regex);

    if(myVariantInfoHG38RegexObj.test(variantIdTerm)){
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

 