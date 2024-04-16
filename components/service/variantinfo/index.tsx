import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/input/term'
import { RegulatoryElementTerm } from '@/components/core/input/term'
import { VariantSet } from '@/components/core/input/set'
import { Table, Cell, Column} from '@/app/components/Table'
import { z } from 'zod'
import { linkeddatahub_icon } from '@/icons'
import { getRegElemPositionData } from '@/components/service/regulatoryElementInfo'
import { downloadBlob } from '@/utils/download'

let caIdRegex = "^(CA|ca)[0-9]";
let rsIdRegex = "^(RS|rs)[0-9]";
let hgvsRegex = "^NM_|XM_[0-9]\.[0-9]:$";
let clinvarRegex = "^[0-9]+$";
let gnomADRegex = "[0-9]*-[0-9]*-[A-Za-z]*-[A-Za-z]*";
let myVariantInfoHG38Regex = "chr[0-9]*:[a-z]\.[0-9]*[A-Za-z]*";

let variantIdResolveErrorMessage = "Unable to resolve the suplied Variant ID to a CAID!";
let alleleRegRespErrorMessage = "Unable to get data from the Allele Registry API, please try again or wait a few minutes before the next attempt!";
let gitDataHubErroMessage = "Unable to get data from Git Data Hub API, please try again or wait a few minutes before the next attempt!";

const AlleleRegistryVariantInfoC = z.object({
  '@id': z.string(),
  entId: z.string(),
  externalRecords: z.any({}).optional()
})
export type AlleleRegistryVariantInfo = z.infer<typeof AlleleRegistryVariantInfoC>

const AlleleRegistryVariantSetInfoC = z.array(
  AlleleRegistryVariantInfoC
);

const AlleleRegistryExternalSourcesInfoC = z.array(
  z.object({
    name: z.string(),
    sources: z.array(z.object({
      '@id':z.string(),
      id: z.string()
    }))
  })
).nullable();
export type AlleleRegistryExternalSourcesInfo = z.infer<typeof AlleleRegistryExternalSourcesInfoC>


const AlleleRegistryExternalSourcesSetInfoC = z.array(
  z.object({
    variantCaId: z.string(),
    externalRecords: AlleleRegistryExternalSourcesInfoC
  })
);

const AlleleSpecificEvidenceInfoC = z.array(
  z.object({
    ldhId: z.string(),
    ldhIri: z.string(),
    sourceDescription: z.string(),
    alleleSpecificityList: z.array(z.object({
      name: z.string(),
      altAlleleQuant: z.any(),
      refAlleleQuant: z.any(),
      sig: z.any()
    }))
  })
);

const xQTL_EvidenceEntContent = z.object({
  GTExIri: z.string(),
  score: z.number(),
  sourceDescription: z.string(),
  eQTL: z.object({
    nes: z.string(),
    sig: z.string()
  }).optional(),
  sQTL: z.object({
    nes: z.string(),
    sig: z.string()
  }).optional(),
  esQTL: z.object({
    nes: z.string(),
    sig: z.string()
  }).optional(),
  type: z.string().optional()
});

const xQTL_EvidenceArray = z.array(
  z.object({
    ldhId: z.string(),
    xQTLEntContent: xQTL_EvidenceEntContent
  })
);

export const GitHubVariantInfoC =  z.object({
  data: z.object({
    entId: z.string(),
    entType: z.string(),
    ld: z.object({
      AlleleSpecificEvidence: z.array(
        z.object({
          ldhId: z.string(),
          ldhIri: z.string(),
          entContent: z.object({
            sourceDescription: z.string(),
            AlleleSpecificity:z.any().optional()
          })
        })
      ),
      xqtlEvidence:  z.array(z.object({
        entId: z.string(),
        ldhId: z.string(),
        ldhIri: z.string(),
        entContent: xQTL_EvidenceEntContent
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
export type GitHubVariantInfo = z.infer<typeof GitHubVariantInfoC>

const HG38SingleGeneAssociationsC = z.object({
  distance_to_feature : z.string().optional(),
  effect : z.string().optional(),
  feature_id : z.string().optional(),
  feature_type : z.string().optional(),
  gene_id : z.string().optional(),
  genename : z.string().optional(),
  hgvs_c : z.string().optional(),
  putative_impact : z.string().optional(),
  transcript_biotype: z.string().optional()
});

const HG38GeneAssociationsC = z.object({
  snpeff: z.object({
    ann: z.array( HG38SingleGeneAssociationsC )
  })  
});
type HG38GeneAssociations = z.infer<typeof HG38GeneAssociationsC>

const HG38GeneAssociationsProcessedForViewC = z.array(
  z.object({
    geneId: z.string(),
    associations: z.array( HG38SingleGeneAssociationsC )
  })
)

const HG38GeneAssociationsSetC = z.array(
  z.object({
    variantCaId: z.string(),
    geneAssociationsHg38: HG38GeneAssociationsProcessedForViewC
  })
)

/*
const GenomeNetworkAlleleObjReponse = z.object({
  '@id':z.string(),
})
const GenomeNetworkAlleleArrayReponse = z.array( GenomeNetworkAlleleObjReponse )
*/

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

async function resolveVarinatCaID(variantIdTerm: string){
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

async function getGitDataHubVariantInfo(variantId: string): Promise<GitHubVariantInfo> {
  const res = await fetch(`https://ldh.genome.network/cfde/ldh/Variant/id/${encodeURIComponent(variantId)}`)
  return await res.json()
}

export async function getAlleleRegistryVariantInfo(variantId: string): Promise<AlleleRegistryVariantInfo> {
  const res = await fetch(`https://reg.genome.network/allele/${encodeURIComponent(variantId)}`);
  return await res.json()
}

export const VariantInfo = MetaNode('VariantInfo')
  .meta({
    label: 'Variant Information',
    description: 'A Variant resolved with reg.clinicalgenome.org',
    icon: [linkeddatahub_icon],
  })
  .codec(AlleleRegistryVariantInfoC)
  .view(variantinfo => (
    <div className="prose max-w-none">
      <a target="_blank" href={`https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=${variantinfo.entId}`}>{variantinfo.entId}</a> (variant)
    </div>
  ))
  .build()

export const VariantInfoFromVariantTerm = MetaNode('VariantInfoFromVariantTerm')
  .meta({
    label: 'Resolve Variant Info from Term',
    description: 'Resolve variant info (Allele registry API) from variant term',
    icon: [linkeddatahub_icon],
  })
  .inputs({ variant: VariantTerm })
  .output(VariantInfo)
  .resolve(async (props) => {
    var varCaId = await resolveVarinatCaID(props.inputs.variant);
    if(varCaId == null || varCaId == ''){
      throw new Error(variantIdResolveErrorMessage);
    }

    console.log("Var CaId: "+varCaId)
    const response = await getAlleleRegistryVariantInfo(varCaId);
    if(response == null){
      throw new Error(alleleRegRespErrorMessage);
    }

    response.entId = varCaId;
    return response;
  })
  .story(props => ``)
  .build()
  
  export const VariantSetInfo = MetaNode('VariantSetInfo')
  .meta({
    label: 'Variant Set Info',
    description: ''
  })
  .codec(AlleleRegistryVariantSetInfoC)
  .view(variantInfoSet => {
    return (
      <Table
        height={500}
        cellRendererDependencies={[variantInfoSet]}
        numRows={variantInfoSet.length}
        enableGhostCells
        enableFocusedCell
      >
        <Column
          name="Variant id (link)"
          cellRenderer={row => <Cell key={row+''}><a target="_blank" href={`https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=${variantInfoSet[row].entId}`}>{variantInfoSet[row].entId}</a></Cell>}
        />
      </Table>
    )
  }).build()

  export const VariantInfoFromVariantSet = MetaNode('VariantInfoFromVariantSet')
  .meta({
    label: `Variant Info From Variant Set`,
    description: "Get Variant Info from Allele Registry for a set of Variant CaId's."
  })
  .inputs({ variantset: VariantSet })
  .output(VariantSetInfo)
  .resolve(async (props) => {
    var varinatSetInpt = props.inputs.variantset.set;
    
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
  }).story(props =>
    `Get Allele Registry data for Variant Set.`
  ).build()

export const GetRegulatoryElementsForThisVariant = MetaNode('GetRegulatoryElementForThisVariant')
  .meta({
    label: 'Identify Variant And Regulatory Element Association',
    description: 'Get Regulatory Elements For This Variant ID',
  })
  .inputs({ variant: VariantTerm })
  .output(RegulatoryElementTerm)
  .resolve(async (props) => {
    var varCaId = await resolveVarinatCaID(props.inputs.variant);
    if(varCaId == null || varCaId == ''){
      throw new Error(variantIdResolveErrorMessage);
    }

    const response = await getGitDataHubVariantInfo(varCaId);
    if(response == null || response.data == null){
      throw new Error(gitDataHubErroMessage);
    }

    if(response.data.ldFor.RegulatoryElement != null && response.data.ldFor.RegulatoryElement.length == 1){
      return response.data.ldFor.RegulatoryElement[0].entId;
    }
    return "N/A";
  })
  .story(props => ``)
  .build()

export const AlleleSpecificEvidencesTable = MetaNode('AlleleSpecificEvidencesTable')
  .meta({
    label: 'Allele Specific Evidences Table',
    description: 'A table of allele specific evidences',
  })
  .codec(AlleleSpecificEvidenceInfoC)
  .view(alleleSpecificEvidence => {
      return (
          <Table
            height={500}
            cellRendererDependencies={[alleleSpecificEvidence]}
            numRows={alleleSpecificEvidence.length}
            enableGhostCells
            downloads={{
              JSON: () => downloadBlob(new Blob([JSON.stringify(alleleSpecificEvidence)], { type: 'application/json;charset=utf-8' }), 'data.json')
            }}
          >
            <Column
              name="LDH Id"
              cellRenderer={row => <Cell key={row+''}>{alleleSpecificEvidence[row].ldhId}</Cell>}
            />
            <Column
              name="Tissue Site or Cell Type"
              cellRenderer={row => <Cell key={row+''}>{alleleSpecificEvidence[row].sourceDescription.replace(/_/g, " ")}</Cell>}
            />
            <Column
              name="LDH Iri"
              cellRenderer={row => <Cell key={row+''}><a target="_blank" href={`${alleleSpecificEvidence[row].ldhIri}`}>evidence link</a></Cell>}
            />
            <Column
              name="Allele Specificity Type"
              cellRenderer={row =>
                <Cell key={row+''}>
                  <table style={{borderCollapse: 'collapse', width:'100%'}}>
                    {alleleSpecificEvidence[row].alleleSpecificityList.map(sources =>
                        <tr><td>{ sources.name }</td></tr>
                    )}
                  </table>
                </Cell>}
            />
            <Column
              name="Allele Specificity Ref. Count"
              cellRenderer={row =>
              <Cell key={row+''}>{
                <table style={{borderCollapse: 'collapse', width:'100%'}}>
                {alleleSpecificEvidence[row].alleleSpecificityList.map(sources =>
                    <tr><td>{ sources.refAlleleQuant }</td></tr>
                  )}
                </table>
              }</Cell>}
            />
            <Column
              name="Allele Specificity Alt. Count"
              cellRenderer={row => <Cell key={row+''}>{
                <table style={{borderCollapse: 'collapse', width:'100%'}}>
                {alleleSpecificEvidence[row].alleleSpecificityList.map(sources =>
                    <tr><td>{ sources.altAlleleQuant }</td></tr>
                  )}
                </table>
              }</Cell>}
            />
            <Column
              name="Adjusted P-Value"
              cellRenderer={row => <Cell key={row+''}>{
                <table style={{borderCollapse: 'collapse', width:'100%'}}>
                {alleleSpecificEvidence[row].alleleSpecificityList.map(sources =>
                    <tr><td>{ sources.sig }</td></tr>
                  )}
                </table>
              }</Cell>}
            />
          </Table>
      )
  })
  .build()

function getAlleleSpecificEvdncFromGitDataHub(alleleSpecificEvidencesList: any){
  var alleleSpecificEvidence: any = [];
  for(let a in alleleSpecificEvidencesList){
    let specificities = alleleSpecificEvidencesList[a].entContent.AlleleSpecificity;
    if(specificities == null){
      continue;
    }
    let specificitieNamesList = Object.getOwnPropertyNames(specificities);
    if(specificitieNamesList.length == 0){
      continue;
    }
    let alleleSpecificityList: any = []

    for(let s in specificitieNamesList){
        let specificitieName = specificitieNamesList[s];
        let specificity = specificities[specificitieName.toString()];

        var specificityObject: any = null;

        if(specificitieName != "HM" && specificitieName != "TF"){
          specificityObject = {
            "name": specificitieName,
            "altAlleleQuant": specificity.altAlleleQuant,
            "refAlleleQuant": specificity.refAlleleQuant,
            "sig": specificity.sig
          }
          alleleSpecificityList.push(specificityObject);
          continue;
        }else{
          let specificitySubName = Object.getOwnPropertyNames(specificity);
          for(let ssN in specificitySubName){
            let s = specificity[specificitySubName[ssN]];
            specificityObject = {
              "name": specificitieName+": "+specificitySubName[ssN],
              "altAlleleQuant": s.altAlleleQuant,
              "refAlleleQuant": s.refAlleleQuant,
              "sig": s.sig
            }
            alleleSpecificityList.push(specificityObject);
            continue;
          }
        }
    }

    let specificEvidence = {
      "ldhId": alleleSpecificEvidencesList[a].ldhId,
      "ldhIri": alleleSpecificEvidencesList[a].ldhIri,
      "sourceDescription": alleleSpecificEvidencesList[a].entContent.sourceDescription,
      "alleleSpecificityList": alleleSpecificityList
    }
    alleleSpecificEvidence.push(specificEvidence)
  }
  return alleleSpecificEvidence;
}

export const GetAlleleSpecificEvidencesForThisVariant = MetaNode('GetAlleleSpecificEvidencesForThisVariant')
  .meta({
    label: 'Resolve Allele Specific Evidences',
    description: 'Get allele specific evidences for this variant',
  })
  .inputs({ variant: VariantTerm })
  .output(AlleleSpecificEvidencesTable)
  .resolve(async (props) => {
    var varCaId = await resolveVarinatCaID(props.inputs.variant);
    if(varCaId == null || varCaId == ''){
      throw new Error(variantIdResolveErrorMessage);
    }

    const response = await getGitDataHubVariantInfo(varCaId);
    if(response == null || response.data == null){
      throw new Error(gitDataHubErroMessage);
    }

    let alleleSpecificEvidencesList = null;
    if(response.data.ld != null &&  response.data.ld.AlleleSpecificEvidence != null){
      alleleSpecificEvidencesList = response.data.ld.AlleleSpecificEvidence;
    }
    return getAlleleSpecificEvdncFromGitDataHub(alleleSpecificEvidencesList);
  })
  .story(props => `Allele specific evidences for the variant${props.inputs ? ` ${props.inputs.variant}` : ''} were resolved.`)
  .build()

export const xQTL_EvidenceDataTable = MetaNode('xQTL_EvidenceDataTable')
  .meta({
    label: 'xQTL Evidence Data Table',
    description: ''
  })
  .codec(GitHubVariantInfoC)
  .view(variantinfo => {
    let xqtlEvidences = variantinfo.data.ld.xqtlEvidence;
    return (
      <Table
        height={500}
        cellRendererDependencies={[xqtlEvidences]}
        numRows={xqtlEvidences.length}
        enableGhostCells
        enableFocusedCell
        downloads={{
          JSON: () => downloadBlob(new Blob([JSON.stringify(xqtlEvidences)], { type: 'application/json;charset=utf-8' }), 'data.json')
        }}
      >
        <Column
          name="LHD Id"
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].ldhId}</Cell>}
        />
        <Column
          name="Evidence link"
          cellRenderer={row => <Cell key={row+''}><a target="_blank" href={`${xqtlEvidences[row].entContent.GTExIri}`}>evidence link</a></Cell>}
        />
        <Column
          name="QTL Type"
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].entContent.type}</Cell>}
        />
        <Column
          name="Normalized Effect Size (nes)"
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].entContent.esQTL?.nes ?? null}</Cell>}
        />
        <Column
          name="P-Value"
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].entContent.esQTL?.sig ?? null}</Cell>}
        />
        <Column
          name="Tissue site"
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].entContent.sourceDescription.replace(/_/g, " ")}</Cell>}
        />
      </Table>
    )
  })
  .build()

function reformatxQTLEvidences(xqtlEvidences: any){
  for(let e in xqtlEvidences){
    let xqtlE_entContent = xqtlEvidences[e].entContent;
    if(xqtlE_entContent.hasOwnProperty('eQTL')){
      xqtlE_entContent.esQTL = xqtlE_entContent.eQTL;
      xqtlE_entContent.type = "eQTL";
      delete xqtlE_entContent.eQTL;
    }else if(xqtlE_entContent.hasOwnProperty('sQTL')){
      xqtlE_entContent.esQTL = xqtlE_entContent.sQTL;
      delete xqtlE_entContent.sQTL;
      xqtlE_entContent.type = "sQTL";
    }
  }
}

function reformatxQTLEvidences2(xqtlEvidences: any){
  let newXQTLArray : any= [];
  for(let e in xqtlEvidences){
    let xqtlE_entContent = xqtlEvidences[e].entContent;
    if(xqtlE_entContent.hasOwnProperty('eQTL')){
      xqtlE_entContent.esQTL = xqtlE_entContent.eQTL;
      xqtlE_entContent.type = "eQTL";
      delete xqtlE_entContent.eQTL;
    }else if(xqtlE_entContent.hasOwnProperty('sQTL')){
      xqtlE_entContent.esQTL = xqtlE_entContent.sQTL;
      delete xqtlE_entContent.sQTL;
      xqtlE_entContent.type = "sQTL";
    }
    let tempObj = {
      'ldhId': xqtlEvidences[e].ldhId,
      'xQTLEntContent': xqtlE_entContent
    }
    newXQTLArray.push(tempObj);
  }
  return newXQTLArray;
}

export const GetxQTL_EvidencesDataForVariantInfo = MetaNode('GetxQTL_EvidencesDataForVariantInfo')
  .meta({
    label: 'Resolve xQTL Evidence Data for Variant Info',
    description: 'Resolve xQTL Evidence Data for Variant Info Data',
  })
  .inputs({ variant: VariantTerm  })
  .output(xQTL_EvidenceDataTable)
  .resolve(async (props) => {
    var varCaId = await resolveVarinatCaID(props.inputs.variant);
    if(varCaId == null || varCaId == ''){
      throw new Error(variantIdResolveErrorMessage);
    }

    console.log("Var CaId: "+varCaId)
    let response = await getGitDataHubVariantInfo(varCaId);
    if(response == null  || response.data == null){
      throw new Error(gitDataHubErroMessage);
    }

    if(response.data.ld != null &&  response.data.ld.xqtlEvidence != null){
      reformatxQTLEvidences(response.data.ld.xqtlEvidence);
    }else{
      throw new Error("Unable to process data, missing values!");
    }
    
    return response;
  })
  .story(props => `xQTL evidence data for the variant${props.inputs ? ` ${props.inputs.variant}` : ''} was resolved.`)
  .build()

export const AlleleRegistryExternalRecordsTable = MetaNode('AlleleRegistryExternalRecordsTable')
  .meta({
    label: 'Allele Registry External Records Table',
    description: ''
  })
  .codec(AlleleRegistryExternalSourcesInfoC)
  .view(AlleleRegistryExternalSourcesList => {
    let sourcesList = AlleleRegistryExternalSourcesList ?? [];

    return (
      <Table
        height={500}
        cellRendererDependencies={[sourcesList]}
        numRows={sourcesList.length}
        enableGhostCells
        enableFocusedCell
        downloads={{
          JSON: () => downloadBlob(new Blob([JSON.stringify(sourcesList)], { type: 'application/json;charset=utf-8' }), 'data.json')
        }}
      >
        <Column
          name="Data Base Name"
          cellRenderer={row => <Cell key={row+''}>{sourcesList[row].name}</Cell>}
        />
        <Column
          name="Variant Id"
          cellRenderer={row =>
          <Cell  key={row+''}>
              <table style={{borderCollapse: 'collapse', width:'100%'}}>
                  {sourcesList[row].sources.map(sources =>
                      <tr><td>{ sources.id }</td></tr>
                  )}
              </table>
          </Cell>}
        />
        <Column
          name="Link"
          cellRenderer={row =>
          <Cell  key={row+''}>
              <table style={{borderCollapse: 'collapse', width:'100%'}}>
                    {sourcesList[row].sources.map(sources =>
                        <tr><td><a target="_blank" href={`${sources['@id']}`}>Link</a></td></tr>
                    )}
              </table>
          </Cell>}
        />
      </Table>
    )
  })
  .build()

function processExternalRecords(variantInfoObj: AlleleRegistryVariantInfo){
  let alleleInfoExternalResources: AlleleRegistryExternalSourcesInfo= [];
  let externalSources = variantInfoObj['externalRecords'];
  for(let er in externalSources){
    if(externalSources[er] != null){

      let extSources = externalSources[er];
      for(let indxEs in extSources){
        var es = extSources[indxEs];
        if(es.id == null && es.rs != null){
          es.id = es.rs.toString();
        }else if(es.id == null && es.preferredName != null){
          //for ClinVarAlleles
          es.id = es.preferredName.toString();
        }else if(es.id == null && es.variationId != null){
          //ClinVarVariations
          es.id = es.variationId.toString();
        }
      }

      let  externalResourcesTemp = {
        name: er.toString(),
        sources: extSources
      };
      alleleInfoExternalResources.push(externalResourcesTemp);
    }
  }
  return alleleInfoExternalResources;
}

export const GetAlleleRegistryExternalRecordsForVariant = MetaNode('GetAlleleRegistryExternalRecordsForVariant')
  .meta({
    label: 'Resolve Allele Registry External Records',
    description: 'Get allele registry external records',
  })
  .inputs({ variant: VariantTerm  })
  .output(AlleleRegistryExternalRecordsTable)
  .resolve(async (props) => {
    var varCaId = await resolveVarinatCaID(props.inputs.variant);
    if(varCaId == null || varCaId == ''){
      throw new Error(variantIdResolveErrorMessage);
    }

    let variantInfoObj: any = await getAlleleRegistryVariantInfo(varCaId);
    if(variantInfoObj == null){
      throw new Error(alleleRegRespErrorMessage);
    }

    let reponse = null;
    if(variantInfoObj['externalRecords'] != null){
      reponse = processExternalRecords(variantInfoObj);
    }
    return reponse;
  })
  .story(props => `External records for the variant${props.inputs ? ` ${props.inputs.variant}` : ''} were resolved.`)
  .build()

export const GeneAssociations_HG38 = MetaNode('GeneAssociations_HG38')
  .meta({
    label: 'Gene Associations HG38',
    description: ''
  })
  .codec(HG38GeneAssociationsProcessedForViewC)
  .view(GeneAssociationsList => {
    return (
      <Table
        height={500}
        cellRendererDependencies={[GeneAssociationsList]}
        numRows={GeneAssociationsList.length}
        enableGhostCells
        enableFocusedCell
        downloads={{
          JSON: () => downloadBlob(new Blob([JSON.stringify(GeneAssociationsList)], { type: 'application/json;charset=utf-8' }), 'data.json')
        }}
      >
        <Column
          name="Gene ID"
          cellRenderer={row => <Cell key={row+''}>{GeneAssociationsList[row].geneId}</Cell>}
        />
        <Column
          name="Variant association to Transcript"
          cellRenderer={row => <Cell key={row+''}>
                <table style={{borderCollapse: 'collapse', width:'100%'}}>
                  {GeneAssociationsList[row].associations.map(associations =>
                      <tr><td>{ associations.effect }</td></tr>
                  )}
                </table>
          </Cell>}
        />
        <Column
          name="Distance To Transcript (bp)"
          cellRenderer={row => <Cell key={row+''}>
                <table style={{borderCollapse: 'collapse', width:'100%'}}>
                  {GeneAssociationsList[row].associations.map(associations =>
                      <tr><td>{ associations.distance_to_feature }</td></tr>
                  )}
                </table>
          </Cell>}
        />
        <Column
          name="Transcript ID"
          cellRenderer={row => <Cell key={row+''}>
                <table style={{borderCollapse: 'collapse', width:'100%'}}>
                  {GeneAssociationsList[row].associations.map(associations =>
                      <tr><td>{ associations.feature_id }</td></tr>
                  )}
                </table>
          </Cell>}
        />
      </Table>
    )
  })
  .build()

  function processHG38externalRecords(externalRecord: any){
    let sourceLinks = externalRecord.sources;
    for(let sl in sourceLinks){
      return sourceLinks[sl]['@id'];
    }
  }

  function processHG38ExternalRecordsResponse(apiResponse: any){
    let geneID: string = "";
    var associatedGeensMap: any = {};
    var associatedGeens = apiResponse.snpeff.ann;
    for(let agIdx in associatedGeens){
      let associatedGeenObj = associatedGeens[agIdx];
      geneID = associatedGeenObj.gene_id+"";

      if(associatedGeenObj.distance_to_feature == null){
        associatedGeenObj['distance_to_feature'] = "within gene";
      }

      if(!associatedGeensMap.hasOwnProperty(geneID) && geneID != ""){
        associatedGeensMap[geneID] = [];
        associatedGeensMap[geneID].push(associatedGeenObj);
      }else if(associatedGeensMap.hasOwnProperty(geneID) && geneID != ""){
        associatedGeensMap[geneID].push(associatedGeenObj);
      }
    }

    let associatedGeensList = [];
    for(let agmInd in associatedGeensMap){
      let associatedGeenFInalData = {
          "geneId": agmInd,
          "associations":associatedGeensMap[agmInd]
      }
      associatedGeensList.push(associatedGeenFInalData)
    }

    return associatedGeensList;
  }

  async function  getHG38externalRecordsDataAPI(hg38ExtSourceLink: string): Promise<HG38GeneAssociations>{
    const req = await fetch(hg38ExtSourceLink, {
      method: 'GET',
      headers: {
        Accept: 'text/plain',
      },
    })
    return await req.json();
  }

  export const GetVarinatToGeneAssociation_HG38 = MetaNode('GetVarinatToGeneAssociation_HG38')
  .meta({
    label: `Get Variant To Gene Association HG38`,
    description: "Get Associated Gene info for a given Variant."
  })
  .inputs({ externalRecords: AlleleRegistryExternalRecordsTable })
  .output(GeneAssociations_HG38)
  .resolve(async (props) => {
    let externalRecords: any = props.inputs.externalRecords;

    let hg38ExtSourceLink = null;
    let response = null;

    for(let er in externalRecords){
      if(externalRecords[er] != null && externalRecords[er].name == "MyVariantInfo_hg38"){
        hg38ExtSourceLink = processHG38externalRecords(externalRecords[er]);
        break;
      }
    }

    if(hg38ExtSourceLink != null){
      response = await getHG38externalRecordsDataAPI(hg38ExtSourceLink);
    }

    return processHG38ExternalRecordsResponse(response);
  }).story(props =>
    `Get Associated Gene info for a given Variant.`
  ).build()

  export const VarinatSetExternalRecordsInfo = MetaNode('VarinatSetExternalRecordsInfo')
  .meta({
    label: 'Variant Set External Records Info',
    description: ''
  })
  .codec(AlleleRegistryExternalSourcesSetInfoC)
  .view( externalRecordsSet => {
    //let externalRecords = externalRecordsSet;
    return ( 
      <Table
      height={500}
      cellRendererDependencies={[externalRecordsSet]}
      numRows={externalRecordsSet.length}
      enableGhostCells
      enableFocusedCell
      downloads={{
        JSON: () => downloadBlob(new Blob([JSON.stringify(externalRecordsSet)], { type: 'application/json;charset=utf-8' }), 'data.json')
      }}
      >
        <Column
          name="Variant CaId"
          cellRenderer={row => <Cell key={row+''}>{externalRecordsSet[row].variantCaId}</Cell>}
        />
        <Column
          name="External Resource Name"
          cellRenderer={row =>
            <Cell key={row+''}>
              <table style={{borderCollapse: 'collapse', width:'100%'}}>
                {externalRecordsSet[row].externalRecords?.map(externalRecord =>
                    <tr><td>{ externalRecord.name}</td></tr>
                )}
              </table>
            </Cell>}
        />
        <Column
          name="Resource ID"
          cellRenderer={row =>
            <Cell key={row+''}>
              <table style={{borderCollapse: 'collapse', width:'100%'}}>
                {externalRecordsSet[row].externalRecords?.map(externalRecord =>
                    <tr><td>{ externalRecord.sources[0].id}</td></tr>
                )}
              </table>
            </Cell>}
        />
        <Column
          name="Resource link"
          cellRenderer={row =>
            <Cell key={row+''}>
              <table style={{borderCollapse: 'collapse', width:'100%'}}>
                {externalRecordsSet[row].externalRecords?.map(externalRecord =>
                    <tr><td><a target="_blank" href={externalRecord.sources[0]['@id']}>Resource link</a></td></tr>
                )}
              </table>
            </Cell>}
          />
      </Table>
    )  
  })
  .build()

  export const GetVarinatSetExternalRecords = MetaNode('GetVarinatSetExternalRecords')
  .meta({
    label: `Get Variant Set External Records`,
    description: "Get External Records for a given Variant Set."
  })
  .inputs({ variantSetInfo: VariantSetInfo })
  .output(VarinatSetExternalRecordsInfo)
  .resolve(async (props) => {
    let variantSetInfo = props.inputs.variantSetInfo;
    let externalRecordsSet = [];

    for(let indx in variantSetInfo){
      let variantInfoObj = variantSetInfo[indx];
      if(variantInfoObj['externalRecords'] != null){
        let exteralRecordsData = processExternalRecords(variantInfoObj);
        let variantExternalRecords = {
          "variantCaId" : variantInfoObj.entId,
          "externalRecords" : exteralRecordsData
        };

        externalRecordsSet.push(variantExternalRecords);
      }
    }
    return externalRecordsSet;
  }).story(props =>
    `Get External Records for a given Variant Set.`
  ).build()

  export const GeneAssociationsSet_HG38 = MetaNode('GeneAssociationsSet_HG38')
  .meta({
    label: 'Gene Associations Set HG38',
    description: ''
  })
  .codec(HG38GeneAssociationsSetC)
  .view( geneAssociationsSet => {
    return ( 
      <Table
      height={500}
      cellRendererDependencies={[geneAssociationsSet]}
      numRows={geneAssociationsSet.length}
      enableGhostCells
      enableFocusedCell
      downloads={{
        JSON: () => downloadBlob(new Blob([JSON.stringify(geneAssociationsSet)], { type: 'application/json;charset=utf-8' }), 'data.json')
      }}
      >
        <Column
          name="Variant CaId"
          cellRenderer={row => <Cell key={row+''}>{ geneAssociationsSet[row].variantCaId }</Cell>}
        />
        <Column
          name="Gene ID"
          cellRenderer={row =>
            <Cell key={row+''}>
              <table style={{borderCollapse: 'collapse', width:'100%'}}>
                {geneAssociationsSet[row].geneAssociationsHg38.map(associationsList =>              
                    <tr style={{borderBottom:'1px solid lightgrey'}}>
                      <td>{ associationsList.geneId }</td>
                      <td>
                        <table style={{borderCollapse: 'collapse', width:'100%'}}>
                            {associationsList.associations.map(association =>
                                <tr><td style={{visibility:'hidden'}}>{association.effect}</td></tr>
                            )}
                        </table>
                      </td>
                    </tr>     
                )}
              </table>
            </Cell>}
        />
        <Column
          name="Effect"
          cellRenderer={row =>
            <Cell key={row+''}>             
                {geneAssociationsSet[row].geneAssociationsHg38.map(associationsList =>
                    <table style={{borderCollapse: 'collapse', borderBottom:'1px solid lightgrey', width:'100%'}}>
                        {associationsList.associations.map(association =>
                            <tr><td>{ association.effect }</td></tr>
                        )}
                    </table>
                )}             
            </Cell>}
        />
        <Column
          name="Feature ID"
          cellRenderer={row =>
            <Cell key={row+''}>             
                {geneAssociationsSet[row].geneAssociationsHg38.map(associationsList =>
                    <table style={{borderCollapse: 'collapse', borderBottom:'1px solid lightgrey', width:'100%'}}>
                        {associationsList.associations.map(association =>
                            <tr><td>{ association.feature_id }</td></tr>
                        )}
                    </table>
                )}             
            </Cell>}
        />
        <Column
          name="Feature Type"
          cellRenderer={row =>
            <Cell key={row+''}>             
                {geneAssociationsSet[row].geneAssociationsHg38.map(associationsList =>
                    <table style={{borderCollapse: 'collapse', borderBottom:'1px solid lightgrey', width:'100%'}}>
                        {associationsList.associations.map(association =>
                            <tr><td>{ association.feature_type }</td></tr>
                        )}
                    </table>
                )}             
            </Cell>}
        />
        <Column
          name="Distance To Feature (bp)"
          cellRenderer={row =>
            <Cell key={row+''}>             
                {geneAssociationsSet[row].geneAssociationsHg38.map(associationsList =>
                    <table style={{borderCollapse: 'collapse', borderBottom:'1px solid lightgrey', width:'100%'}}>
                        {associationsList.associations.map(association =>
                            <tr><td>{ association.distance_to_feature }</td></tr>
                        )}
                    </table>
                )}             
            </Cell>}
        />
        <Column
          name="HGVS"
          cellRenderer={row =>
            <Cell key={row+''}>             
                {geneAssociationsSet[row].geneAssociationsHg38.map(associationsList =>
                    <table style={{borderCollapse: 'collapse', borderBottom:'1px solid lightgrey', width:'100%'}}>
                        {associationsList.associations.map(association =>
                            <tr><td>{ association.hgvs_c }</td></tr>
                        )}
                    </table>
                )}             
            </Cell>}
        />
        <Column
          name="Putative Impact"
          cellRenderer={row =>
            <Cell key={row+''}>             
                {geneAssociationsSet[row].geneAssociationsHg38.map(associationsList =>
                    <table style={{borderCollapse: 'collapse', borderBottom:'1px solid lightgrey', width:'100%'}}>
                        {associationsList.associations.map(association =>
                            <tr><td>{ association.putative_impact }</td></tr>
                        )}
                    </table>
                )}             
            </Cell>}
        />
        <Column
          name="Transcript Biotype"
          cellRenderer={row =>
            <Cell key={row+''}>             
                {geneAssociationsSet[row].geneAssociationsHg38.map(associationsList =>
                    <table style={{borderCollapse: 'collapse', borderBottom:'1px solid lightgrey', width:'100%'}}>
                        {associationsList.associations.map(association =>
                            <tr><td>{ association.transcript_biotype }</td></tr>
                        )}
                    </table>
                )}             
            </Cell>}
        />
      </Table>
    )  
  })
  .build()

  export const GetVarinatSetToGeneAssociation_HG38 = MetaNode('GetVarinatSetToGeneAssociation_HG38')
  .meta({
    label: `Identify Variant And Gene Associations`,
    description: "Get Associated Gene info for a given set of Variant External records."
  })
  .inputs({ variantSetExternalRecordsInfo: VarinatSetExternalRecordsInfo })
  .output(GeneAssociationsSet_HG38)
  .resolve(async (props) => {
    let variantExternalRecordsSetInfo = props.inputs.variantSetExternalRecordsInfo;
    let varinatGeneAssociationsSet = [];

    for(let vERIdx in variantExternalRecordsSetInfo){
      let variantExternalRecords = variantExternalRecordsSetInfo[vERIdx];
      let externalRecords: any = variantExternalRecords.externalRecords;
      for(let erIdx in externalRecords){
        if(externalRecords[erIdx] != null && externalRecords[erIdx].name == "MyVariantInfo_hg38"){
          let hg38ExtSourceLink = processHG38externalRecords(externalRecords[erIdx]);
          let response = await getHG38externalRecordsDataAPI(hg38ExtSourceLink);
          let associatedGeensList = processHG38ExternalRecordsResponse(response);

          let varinatGeneAssociation = {
            "variantCaId": variantExternalRecords.variantCaId,
            "geneAssociationsHg38": associatedGeensList
          }
          varinatGeneAssociationsSet.push(varinatGeneAssociation);
        }
      }   
    }

    return varinatGeneAssociationsSet;
  }).story(props =>
    `Get Associated Gene info for a given set of Variant External records.`
  ).build()

  export const REforVariantSetInfo = MetaNode('REforVariantSetInfo')
  .meta({
    label: 'Regulatory Elements For Variant Set Info',
    description: ''
  })
  .codec(
    z.array(z.object({
      caid: z.string(),
      reEntId: z.string(),
      rePosition: z.string()
    }))
  )
  .view( varAndReIdArray => {
    return ( 
      <Table
      height={500}
      cellRendererDependencies={[varAndReIdArray]}
      numRows={varAndReIdArray.length}
      enableGhostCells
      enableFocusedCell
      downloads={{
        JSON: () => downloadBlob(new Blob([JSON.stringify(varAndReIdArray)], { type: 'application/json;charset=utf-8' }), 'data.json')
      }}
      >
        <Column
          name="Variant CaId"
          cellRenderer={row => <Cell key={row+''}>{ varAndReIdArray[row].caid }</Cell>}
        />
        <Column
          name="Regulatory Element EntId"
          cellRenderer={row => <Cell key={row+''}>{ varAndReIdArray[row].reEntId }</Cell>}
        />
        <Column
          name="Regulatory Element Position"
          cellRenderer={row => <Cell key={row+''}>{ varAndReIdArray[row].rePosition }</Cell>}
        />
      </Table>  
    )
  })
  .build()

  export const GetRegulatoryElementsForVariantSet = MetaNode('GetRegulatoryElementsForVariantSet')
  .meta({
    label: `Regulatory Elements For Variant Set`,
    description: "Get Regulatory Elements for the Variant Set"
  })
  .inputs({ variantset: VariantSet })
  .output(REforVariantSetInfo)
  .resolve(async (props) => {
    var varinatSetInpt = props.inputs.variantset.set;
    
    let varAndRegulatoryElem = [];

    for(let indx in varinatSetInpt){
      let varCaID = varinatSetInpt[indx];
      const response = await getGitDataHubVariantInfo(varCaID);
      if(response == null || response.data == null){
        continue;
      }

      if(response.data.ldFor.RegulatoryElement != null && response.data.ldFor.RegulatoryElement.length == 1){
        let reEntId = response.data.ldFor.RegulatoryElement[0].entId;
        let tempObj = {
          'caid':varCaID,
          'reEntId': reEntId,
          'rePosition':""
        }

        const rePositionData = await getRegElemPositionData(reEntId);
        if(rePositionData != null && rePositionData.data.cCREQuery[0].coordinates != null){
          tempObj.rePosition = rePositionData.data.cCREQuery[0].coordinates.chromosome+": "+rePositionData.data.cCREQuery[0].coordinates.start+"-"+rePositionData.data.cCREQuery[0].coordinates.end;
        }
        varAndRegulatoryElem.push(tempObj);
      }
    }

    if(varAndRegulatoryElem.length == 0){
      throw new Error(gitDataHubErroMessage);
    }

    return varAndRegulatoryElem;
  }).story(props =>
    `Get Regulatory Elements id and posotion for Variant Set.`
  ).build()

  export const AlleleSpecificEvidenceForVariantSet = MetaNode('AlleleSpecificEvidenceForVariantSet')
  .meta({
    label: 'Allele Specific Evidence For Variant Set',
    description: ''
  })
  .codec(
    z.array(z.object({
      caid: z.string(),
      alleleSpecificEvidence: AlleleSpecificEvidenceInfoC
    }))
  )
  .view( alleleEvidncForVarSetArray => {
    return ( 
      <Table
      height={500}
      cellRendererDependencies={[alleleEvidncForVarSetArray]}
      numRows={alleleEvidncForVarSetArray.length}
      enableGhostCells
      enableFocusedCell
      downloads={{
        JSON: () => downloadBlob(new Blob([JSON.stringify(alleleEvidncForVarSetArray)], { type: 'application/json;charset=utf-8' }), 'data.json')
      }}
      >
        <Column
          name="Variant CaID"
          cellRenderer={row => <Cell key={row+''}>{alleleEvidncForVarSetArray[row].caid}</Cell>}
        />
        <Column
          name="LDH Id"
          cellRenderer={row =>
          <Cell key={row+''}>{
            <table style={{borderCollapse: 'collapse', width:'100%'}}>
            {alleleEvidncForVarSetArray[row].alleleSpecificEvidence.map(sources =>
                <tr><td>{ sources.ldhId }</td></tr>
              )}
            </table>
          }</Cell>}
        />
        <Column
          name="Tissue Site or Cell Type"
          cellRenderer={row =>
          <Cell key={row+''}>{
            <table style={{borderCollapse: 'collapse', width:'100%'}}>
            {alleleEvidncForVarSetArray[row].alleleSpecificEvidence.map(sources =>
                <tr><td>{ sources.sourceDescription.replace(/_/g, " ") }</td></tr>
              )}
            </table>
          }</Cell>}
        />
        <Column
          name="LDH Iri"
          cellRenderer={row =>
          <Cell key={row+''}>{
            <table style={{borderCollapse: 'collapse', width:'100%'}}>
            {alleleEvidncForVarSetArray[row].alleleSpecificEvidence.map(sources =>
                <tr><td><a target="_blank" href={`${sources.ldhIri}`}>evidence link</a></td></tr>
              )}
            </table>
          }</Cell>}
        />
        <Column
          name="Allele Specificity Type"
          cellRenderer={row =>
          <Cell key={row+''}>{
            <table style={{borderCollapse: 'collapse', width:'100%'}}>
            {alleleEvidncForVarSetArray[row].alleleSpecificEvidence.map(sources =>
                <tr><td>
                  <table style={{borderCollapse: 'collapse', width:'100%'}}>
                      {  sources.alleleSpecificityList.map(sources =>  
                        <tr><td>{ sources.name }</td></tr>  
                      )}
                  </table>
                </td></tr>
              )}
            </table>
          }</Cell>}
        />
        <Column
          name="Allele Specificity Ref. Count"
          cellRenderer={row =>
          <Cell key={row+''}>{
            <table style={{borderCollapse: 'collapse', width:'100%'}}>
            {alleleEvidncForVarSetArray[row].alleleSpecificEvidence.map(sources =>
                <tr><td>
                  <table style={{borderCollapse: 'collapse', width:'100%'}}>
                      {  sources.alleleSpecificityList.map(sources =>  
                        <tr><td>{ sources.refAlleleQuant }</td></tr>  
                      )}
                  </table>
                </td></tr>
              )}
            </table>
          }</Cell>}
        />
        <Column
          name="Allele Specificity Alt. Count"
          cellRenderer={row =>
          <Cell key={row+''}>{
            <table style={{borderCollapse: 'collapse', width:'100%'}}>
            {alleleEvidncForVarSetArray[row].alleleSpecificEvidence.map(sources =>
                <tr><td>
                  <table style={{borderCollapse: 'collapse', width:'100%'}}>
                      {  sources.alleleSpecificityList.map(sources =>  
                        <tr><td>{ sources.altAlleleQuant }</td></tr>  
                      )}
                  </table>
                </td></tr>
              )}
            </table>
          }</Cell>}
        />
        <Column
          name="Adjusted P-Value"
          cellRenderer={row =>
          <Cell key={row+''}>{
            <table style={{borderCollapse: 'collapse', width:'100%'}}>
            {alleleEvidncForVarSetArray[row].alleleSpecificEvidence.map(sources =>
                <tr><td>
                  <table style={{borderCollapse: 'collapse', width:'100%'}}>
                      {  sources.alleleSpecificityList.map(sources =>  
                        <tr><td>{ sources.sig }</td></tr>  
                      )}
                  </table>
                </td></tr>
              )}
            </table>
          }</Cell>}
        />
      </Table>
    )
  })
  .build()

  export const GetVarinatSetAlleleSpecificEvidence = MetaNode('GetVarinatSetAlleleSpecificEvidence')
  .meta({
    label: `Get Variant Set Allele Specific Evidence`,
    description: "Get Allele Specific Evidence form Variant Set."
  })
  .inputs({ variantSetInfo: VariantSetInfo })
  .output(AlleleSpecificEvidenceForVariantSet)
  .resolve(async (props) => {
    let variantSetInfo = props.inputs.variantSetInfo;

    let varinatSetAlleleSpecificEvdnc = [];
    for(let indx in variantSetInfo){
      let variantInfoObj = variantSetInfo[indx];
      const response = await getGitDataHubVariantInfo(variantInfoObj.entId);
      if(response == null || response.data == null || response.data.ld == null || response.data.ld.AlleleSpecificEvidence == null){
        continue;
      }

      let tempObj = {
        'caid': variantInfoObj.entId,
        'alleleSpecificEvidence' : [] 
      };

      let evidenceReponse = getAlleleSpecificEvdncFromGitDataHub(response.data.ld.AlleleSpecificEvidence);
      if(evidenceReponse == null || evidenceReponse.length == 0){
        continue;
      }
      tempObj.alleleSpecificEvidence = evidenceReponse;
      
      varinatSetAlleleSpecificEvdnc.push(tempObj);
    }

    if(varinatSetAlleleSpecificEvdnc.length == 0){
      throw new Error(gitDataHubErroMessage);
    }

    return varinatSetAlleleSpecificEvdnc;
  }).story(props =>
    `Get Allele Specific Evidence form Variant Set.`
  ).build()


  export const xQTL_EvidenceFroVariantSet = MetaNode('xQTL_EvidenceFroVariantSet')
  .meta({
    label: 'xQTL Evidence For Variant Set',
    description: ''
  })
  .codec(
    z.array(z.object({
      caid: z.string(),
      xQTL_Evidence: xQTL_EvidenceArray
    }))
  )
  .view( xQTLEvdVariantSet => {
    return (
      <Table
        height={500}
        cellRendererDependencies={[xQTLEvdVariantSet]}
        numRows={xQTLEvdVariantSet.length}
        enableGhostCells
        enableFocusedCell
        downloads={{
          JSON: () => downloadBlob(new Blob([JSON.stringify(xQTLEvdVariantSet)], { type: 'application/json;charset=utf-8' }), 'data.json')
        }}
      >
        <Column
          name="Variant CaID"
          cellRenderer={row => <Cell key={row+''}>{xQTLEvdVariantSet[row].caid}</Cell>}
        />
        <Column
          name="LDH Id"
          cellRenderer={row =>
          <Cell key={row+''}>{
            <table style={{borderCollapse: 'collapse', width:'100%'}}>
              {xQTLEvdVariantSet[row].xQTL_Evidence.map(sources =>
                <tr><td>{ sources.ldhId }</td></tr>
              )}
            </table>
          }</Cell>}
        />
        <Column
          name="Evidence link"
          cellRenderer={row =>
          <Cell key={row+''}>{
            <table style={{borderCollapse: 'collapse', width:'100%'}}>
              {xQTLEvdVariantSet[row].xQTL_Evidence.map(sources =>
                <tr><td>{ sources.ldhId }</td></tr>
              )}
            </table>
          }</Cell>}
        />
        <Column
          name="QTL Type"
          cellRenderer={row =>
          <Cell key={row+''}>{
            <table style={{borderCollapse: 'collapse', width:'100%'}}>
              {xQTLEvdVariantSet[row].xQTL_Evidence.map(sources =>
                <tr><td>{ sources.xQTLEntContent.type }</td></tr>
              )}
            </table>
          }</Cell>}
        />
        <Column
          name="Normalized Effect Size (nes)"
          cellRenderer={row =>
          <Cell key={row+''}>{
            <table style={{borderCollapse: 'collapse', width:'100%'}}>
              {xQTLEvdVariantSet[row].xQTL_Evidence.map(sources =>
                <tr><td>{ sources.xQTLEntContent.esQTL?.nes ?? null }</td></tr>
              )}
            </table>
          }</Cell>}
        />
        <Column
          name="P-Value"
          cellRenderer={row =>
          <Cell key={row+''}>{
            <table style={{borderCollapse: 'collapse', width:'100%'}}>
              {xQTLEvdVariantSet[row].xQTL_Evidence.map(sources =>
                <tr><td>{ sources.xQTLEntContent.esQTL?.sig ?? null }</td></tr>
              )}
            </table>
          }</Cell>}
        />
        <Column
          name="Tissue site"
          cellRenderer={row =>
          <Cell key={row+''}>{
            <table style={{borderCollapse: 'collapse', width:'100%'}}>
              {xQTLEvdVariantSet[row].xQTL_Evidence.map(sources =>
                <tr><td>{ sources.xQTLEntContent.sourceDescription.replace(/_/g, " ") }</td></tr>
              )}
            </table>
          }</Cell>}
        />
      </Table>
    )
  }).build()

  export const GetVarinatSetXQTLEvidence = MetaNode('GetVarinatSetXQTLEvidence')
  .meta({
    label: `Get Variant Set xQTL Evidence`,
    description: "Get xQTL Evidence form Variant Set."
  })
  .inputs({ variantSetInfo: VariantSetInfo })
  .output(xQTL_EvidenceFroVariantSet)
  .resolve(async (props) => {
    let variantSetInfo = props.inputs.variantSetInfo;

    let varinatSetXQTLEvidnc = [];
    for(let indx in variantSetInfo){
      let variantInfoObj = variantSetInfo[indx];
      const response = await getGitDataHubVariantInfo(variantInfoObj.entId);
      if(response == null || response.data == null || response.data.ld == null || response.data.ld.xqtlEvidence == null){
        continue;
      }

      let tempObj = {
        'caid': variantInfoObj.entId,
        'xQTL_Evidence': []
      };
    
      let evidenceReponse = reformatxQTLEvidences2(response.data.ld.xqtlEvidence);
      if(evidenceReponse == null || evidenceReponse.length == 0){
        continue;
      }
      tempObj.xQTL_Evidence = evidenceReponse;     

      varinatSetXQTLEvidnc.push(tempObj);
    }

    if(varinatSetXQTLEvidnc.length == 0){
      throw new Error(gitDataHubErroMessage);
    }

    return varinatSetXQTLEvidnc;
  }).story(props =>
    "Get xQTL Evidence form Variant Set."
  ).build()

