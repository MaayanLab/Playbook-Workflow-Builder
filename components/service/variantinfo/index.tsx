import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/input/term'
import { RegulatoryElementTerm } from '@/components/core/input/term'
import { VariantSet } from '@/components/core/input/set'
import { Table, Cell, Column} from '@/app/components/Table'
import { z } from 'zod'
import { linkeddatahub_icon } from '@/icons'
import { downloadBlob } from '@/utils/download'
//import internal from 'stream'

let caIdRegex = "^(CA|ca)[0-9]";
let rsIdRegex = "^(RS|rs)[0-9]";
let hgvsRegex = "^NM_|XM_[0-9]\.[0-9]:$";
let clinvarRegex = "^[0-9]+$";
let gnomADRegex = "[0-9]*-[0-9]*-[A-Za-z]*-[A-Za-z]*";
let myVariantInfoHG38Regex = "chr[0-9]*:[a-z]\.[0-9]*[A-Za-z]*";

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
    variantCId: z.string(),
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
      xqtlEvidence: z.array(
        z.object({
          entId: z.string(),
          ldhId: z.string(),
          ldhIri: z.string(),
          entContent: z.object({
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
          })
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
    variantCId: z.string(),
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
  //console.log(JSON.stringify(jsonObj));

  let alleleRegistrylink = null;
  if(Array.isArray(jsonObj)){
    //console.log("size: "+jsonObj.length);
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
    console.log("dbSNP (RsId)");
    return await getGenomeNetworkAlleles_dbSNP(variantIdTerm.slice(2));
  }else if(hgvsRegexObj.test(variantIdTerm)){
    console.log("HGVS");
    return await getGenomeNetworkAlleles_HGVS(variantIdTerm);
  }else if(clinvarRegexObj.test(variantIdTerm)){
    console.log("ClinVar");
    return await getGenomeNetworkAlleles_ClinVar(variantIdTerm);
  }else if(gnomADRegexObj.test(variantIdTerm)){
    console.log("gnomAD");
    return await getGenomeNetworkAlleles_gnomAD(variantIdTerm);
  }else if(myVariantInfoHG38RegexObj.test(variantIdTerm)){
    console.log("myGenomeHG38");
    return await getGenomeNetworkAlleles_myVariantInfoHG38(variantIdTerm);
  }else{
    return variantIdTerm;
  }
}

async function getGitDataHubVariantInfo(variantId: string): Promise<GitHubVariantInfo> {
  const res = await fetch(`https://genboree.org/cfde-gene-dev/Variant/id/${encodeURIComponent(variantId)}`)
  return await res.json()
}

export async function getAlleleRegistryVariantInfo(variantId: string): Promise<AlleleRegistryVariantInfo> {
  const res = await fetch(`https://reg.genome.network/allele/${encodeURIComponent(variantId)}`);
  return await res.json()
}

export const VariantInfo = MetaNode('VariantInfo')
  .meta({
    label: 'Variant Information',
    description: 'A Variant resolved with MyVariantInfo',
    icon: [linkeddatahub_icon],
  })
  .codec(AlleleRegistryVariantInfoC)
  .view(variantinfo => (
    <div>
      <a target="_blank" href={`https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=${variantinfo.entId}`}>{variantinfo.entId}</a> (variant)
    </div>
  ))
  .build()

export const VariantInfoFromVariantTerm = MetaNode('VariantInfoFromVariantTerm')
  .meta({
    label: 'Resolve Variant Info from Term',
    description: 'Resolve variant info (Allele registry API) from variant term with MyVariantInfo',
    icon: [linkeddatahub_icon],
  })
  .inputs({ variant: VariantTerm })
  .output(VariantInfo)
  .resolve(async (props) => {
    var varCaId = await resolveVarinatCaID(props.inputs.variant);
    console.log("Var CaId: "+varCaId)
    const response = await getAlleleRegistryVariantInfo(varCaId);
    response.entId = varCaId;
    return response;
  })
  .build()
  
  export const VariantSetInfo = MetaNode('VariantSetInfo')
  .meta({
    label: 'VariantSetInfo',
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
          name="Varinat id (link)"
          cellRenderer={row => <Cell key={row+''}><a target="_blank" href={`https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=${variantInfoSet[row].entId}`}>{variantInfoSet[row].entId}</a></Cell>}
        />
      </Table>
    )
  }).build()

  export const VariantInfoFromVariantSet = MetaNode('VariantInfoFromVariantSet')
  .meta({
    label: `VariantInfoFromVariantSet`,
    description: "Get Variant Info from Allele Registry for a set of Varinat CId's."
  })
  .inputs({ variantset: VariantSet })
  .output(VariantSetInfo)
  .resolve(async (props) => {
    var varinatSetInpt = props.inputs.variantset.set;
    
    var varInfoSet = [];
    var varAlleleRegResponse = null;
    for(let indx in varinatSetInpt){
      varAlleleRegResponse = await getAlleleRegistryVariantInfo(varinatSetInpt[indx].trim());
      varAlleleRegResponse.entId = varinatSetInpt[indx].trim();
      varInfoSet.push(varAlleleRegResponse);
    }
    return varInfoSet;
  }).story(props =>
    `Get Allele Registry data for Varinat Set.`
  ).build()

export const GetRegulatoryElementsForThisVariant = MetaNode('GetRegulatoryElementForThisVariant')
  .meta({
    label: 'Resolve Reg. Element from Var. Info',
    description: 'GetRegulatoryElementsForThisVariant',
  })
  .inputs({ variant: VariantTerm })
  .output(RegulatoryElementTerm)
  .resolve(async (props) => {
    var varCaId = await resolveVarinatCaID(props.inputs.variant);
    console.log("Var CaId: "+varCaId)
    const response = await getGitDataHubVariantInfo(varCaId);
    if(response.data.ldFor.RegulatoryElement != null && response.data.ldFor.RegulatoryElement.length == 1){
      return response.data.ldFor.RegulatoryElement[0].entId;
    }
    return "N/A";
  })
  .build()

export const AlleleSpecificEvidencesTable = MetaNode('AlleleSpecificEvidencesTable')
  .meta({
    label: 'AlleleSpecificEvidencesTable',
    description: ''
  })
  .codec(AlleleSpecificEvidenceInfoC)
  .view(alleleSpecificEvidence => {
      return (
          <Table
            height={500}
            cellRendererDependencies={[alleleSpecificEvidence]}
            numRows={alleleSpecificEvidence.length}
            enableGhostCells
            enableFocusedCell
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

export const GetAlleleSpecificEvidencesForThisVariant = MetaNode('GetAlleleSpecificEvidencesForThisVariant')
  .meta({
    label: 'Resolve Allele Specific Evidences from Var. Info',
    description: 'GetAlleleSpecificEvidencesForThisVariant',
  })
  .inputs({ variant: VariantTerm })
  .output(AlleleSpecificEvidencesTable)
  .resolve(async (props) => {
    var varCaId = await resolveVarinatCaID(props.inputs.variant);
    console.log("Var CaId: "+varCaId)
    const response = await getGitDataHubVariantInfo(varCaId);

    var alleleSpecificEvidence: any = [];

    let alleleSpecificEvidencesList = response.data.ld.AlleleSpecificEvidence;
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
  })
  .build()

export const xQTL_EvidenceDataTable = MetaNode('xQTL_EvidenceDataTable')
  .meta({
    label: 'xQTL_EvidenceDataTable',
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

export const GetxQTL_EvidencesDataForVariantInfo = MetaNode('GetxQTL_EvidencesDataForVariantInfo')
  .meta({
    label: 'Resolve xQTL Evidence Data for Variant Info',
    description: 'Resolve xQTL Evidence Data for Variant Info Data',
  })
  .inputs({ variant: VariantTerm })
  .output(xQTL_EvidenceDataTable)
  .resolve(async (props) => {
    var varCaId = await resolveVarinatCaID(props.inputs.variant);
    console.log("Var CaId: "+varCaId)
    let response = await getGitDataHubVariantInfo(varCaId);

    let xqtlEvidences = response.data.ld.xqtlEvidence;
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
    return response;
  })
  .build()

export const AlleleRegistryExternalRecordsTable = MetaNode('AlleleRegistryExternalRecordsTable')
  .meta({
    label: 'AlleleRegistryExternalRecordsTable',
    description: ''
  })
  .codec(AlleleRegistryExternalSourcesInfoC)
  .view(AlleleRegistryExternalSourcesList => {
    let sourcesList = AlleleRegistryExternalSourcesList;

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
    label: 'Resolve Allele Registry External Records for Variant',
    description: 'GetAlleleRegistryExternalRecordsForVariant',
  })
  .inputs({ variant: VariantTerm  })
  .output(AlleleRegistryExternalRecordsTable)
  .resolve(async (props) => {
    var varCaId = await resolveVarinatCaID(props.inputs.variant);
    console.log("Var CaId: "+varCaId)
    let variantInfoObj: any = await getAlleleRegistryVariantInfo(varCaId);
    let reponse = null;
    if(variantInfoObj['externalRecords'] != null){
      reponse = processExternalRecords(variantInfoObj);
    }
    return reponse;
}).build()

  export const GeneAssociations_HG38 = MetaNode('GeneAssociations_HG38')
  .meta({
    label: 'GeneAssociations_HG38',
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
          name="Gene ID/Name"
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
        associatedGeenObj['distance_to_feature'] = "/";
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
    label: `GetVarinatToGeneAssociation_HG38`,
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
    label: 'VarinatSetExternalRecordsInfo',
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
          name="Varinat CId"
          cellRenderer={row => <Cell key={row+''}>{externalRecordsSet[row].variantCId}</Cell>}
        />
        <Column
          name="External Resource Name"
          cellRenderer={row =>
            <Cell key={row+''}>
              <table style={{borderCollapse: 'collapse', width:'100%'}}>
                {externalRecordsSet[row].externalRecords.map(externalRecord =>
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
                {externalRecordsSet[row].externalRecords.map(externalRecord =>
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
                {externalRecordsSet[row].externalRecords.map(externalRecord =>
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
    label: `GetVarinatSetExternalRecords`,
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
          "variantCId" : variantInfoObj.entId,
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
    label: 'GeneAssociationsSet_HG38',
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
          name="Variant CId"
          cellRenderer={row => <Cell key={row+''}>{ geneAssociationsSet[row].variantCId }</Cell>}
        />
        <Column
          name="Gene ID/Name"
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
          name="Feature_Type"
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
          name="HGVS_C"
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
    label: `GetVarinatSetToGeneAssociation_HG38`,
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
            "variantCId": variantExternalRecords.variantCId,
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

