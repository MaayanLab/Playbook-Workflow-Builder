import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/input/term'
import { RegulatoryElementTerm } from '@/components/core/input/term'
import { VariantSet } from '@/components/core/input/set'
import { Table, Cell, Column} from '@/app/components/Table'
import { z } from 'zod'
import { linkeddatahub_icon } from '@/icons'

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
);

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

async function getGitDataHubVariantInfo(variantId: string): Promise<GitHubVariantInfo> {
  const res = await fetch(`https://genboree.org/cfde-gene-dev/Variant/id/${encodeURIComponent(variantId)}`)
  return await res.json()
}

export async function getAlleleRegistryVariantInfo(variantId: string): Promise<AlleleRegistryVariantInfo> {
  const res = await fetch(`https://reg.genome.network/allele/${encodeURIComponent(variantId)}`);
  return await res.json()
}

export const VariantInfoFromVariantTerm = MetaNode('VariantInfoFromVariantTerm')
  .meta({
    label: 'Resolve Variant Info from Term',
    description: 'Resolve variant info (Allele registry API) from variant term with MyVariantInfo',
    icon: [linkeddatahub_icon],
  })
  .inputs({ variant: VariantTerm })
  .output(VariantInfo)
  .resolve(async (props) => {
    const response = await getAlleleRegistryVariantInfo(props.inputs.variant);
    response.entId = props.inputs.variant;
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
          name="Varinat id"
          cellRenderer={row => <Cell key={row+''}>{variantInfoSet[row].entId}</Cell>}
        />
        <Column
          name="Variant link"
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
  .inputs({ variantInfo: VariantInfo  })
  .output(RegulatoryElementTerm)
  .resolve(async (props) => {
    const reponse = await getGitDataHubVariantInfo(props.inputs.variantInfo.entId);
    if(reponse.data.ldFor.RegulatoryElement != null && reponse.data.ldFor.RegulatoryElement.length == 1){
      return reponse.data.ldFor.RegulatoryElement[0].entId;
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
          >
            <Column
              name="ldhId"
              cellRenderer={row => <Cell key={row+''}>{alleleSpecificEvidence[row].ldhId}</Cell>}
            />
            <Column
              name="Tissue Site od Cell Type"
              cellRenderer={row => <Cell key={row+''}>{alleleSpecificEvidence[row].sourceDescription.replace(/_/g, " ")}</Cell>}
            />
            <Column
              name="ldhIri"
              cellRenderer={row => <Cell key={row+''}><a target="_blank" href={`${alleleSpecificEvidence[row].ldhIri}`}>evidence link</a></Cell>}
            />
            <Column
              name="Allele Specif. Name"
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
              name="Allele Specif. Ref. Quant"
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
              name="Allele Specif. Alt. Quant"
              cellRenderer={row => <Cell key={row+''}>{
                <table style={{borderCollapse: 'collapse', width:'100%'}}>
                {alleleSpecificEvidence[row].alleleSpecificityList.map(sources =>
                    <tr><td>{ sources.altAlleleQuant }</td></tr>
                  )}
                </table>
              }</Cell>}
            />
            <Column
              name="Allele Specif. sig"
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
  .inputs({ variantInfo: VariantInfo })
  .output(AlleleSpecificEvidencesTable)
  .resolve(async (props) => {
    let response = await getGitDataHubVariantInfo(props.inputs.variantInfo.entId);

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

    console.log(JSON.stringify(alleleSpecificEvidence));

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
      >
        <Column
          name="ldhId"
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].ldhId}</Cell>}
        />
        <Column
          name="Evidence link"
          cellRenderer={row => <Cell key={row+''}><a target="_blank" href={`${xqtlEvidences[row].entContent.GTExIri}`}>evidence link</a></Cell>}
        />
        <Column
          name="Type"
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].entContent.type}</Cell>}
        />
        <Column
          name="Normalized Effect Size (nes)"
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].entContent.esQTL?.nes ?? null}</Cell>}
        />
        <Column
          name="p-value"
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
  .inputs({ variantInfo: VariantInfo  })
  .output(xQTL_EvidenceDataTable)
  .resolve(async (props) => {
    const reponse = await getGitDataHubVariantInfo(props.inputs.variantInfo.entId);
    let xqtlEvidences = reponse.data.ld.xqtlEvidence;
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
    return reponse;
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
        cellRendererDependencies={[AlleleRegistryExternalSourcesList]}
        numRows={AlleleRegistryExternalSourcesList.length}
        enableGhostCells
        enableFocusedCell
      >
        <Column
          name="External Source Name"
          cellRenderer={row => <Cell key={row+''}>{AlleleRegistryExternalSourcesList[row].name}</Cell>}
        />
        <Column
          name="Source Id"
          cellRenderer={row =>
          <Cell  key={row+''}>
              <table style={{borderCollapse: 'collapse', width:'100%'}}>
                  {AlleleRegistryExternalSourcesList[row].sources.map(sources =>
                      <tr><td>{ sources.id }</td></tr>
                  )}
              </table>
          </Cell>}
        />
        <Column
          name="Source Link"
          cellRenderer={row =>
          <Cell  key={row+''}>
              <table style={{borderCollapse: 'collapse', width:'100%'}}>
                    {AlleleRegistryExternalSourcesList[row].sources.map(sources =>
                        <tr><td><a target="_blank" href={`${sources['@id']}`}>Resource link</a></td></tr>
                    )}
              </table>
          </Cell>}
        />
      </Table>
    )
  })
  .build()

  function processExternalRecords(variantInfoObj: AlleleRegistryVariantInfo){
    let alleleInfoExternalResources = [];
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
  .inputs({ variantInfo: VariantInfo  })
  .output(AlleleRegistryExternalRecordsTable)
  .resolve(async (props) => {
    let variantInfoObj: any = props.inputs.variantInfo;
    if(variantInfoObj['externalRecords'] != null){
      return processExternalRecords(variantInfoObj);
    }
    return null;
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
      >
        <Column
          name="Gene ID/Name"
          cellRenderer={row => <Cell key={row+''}>{GeneAssociationsList[row].geneId}</Cell>}
        />
        <Column
          name="Effect"
          cellRenderer={row => <Cell key={row+''}>
                <table style={{borderCollapse: 'collapse', width:'100%'}}>
                  {GeneAssociationsList[row].associations.map(associations =>
                      <tr><td>{ associations.effect }</td></tr>
                  )}
                </table>
          </Cell>}
        />
        <Column
          name="Distance To Feature"
          cellRenderer={row => <Cell key={row+''}>
                <table style={{borderCollapse: 'collapse', width:'100%'}}>
                  {GeneAssociationsList[row].associations.map(associations =>
                      <tr><td>{ associations.distance_to_feature }</td></tr>
                  )}
                </table>
          </Cell>}
        />
        <Column
          name="Feature ID"
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
    let externalRecords = props.inputs.externalRecords;

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
      >
        <Column
          name="Varinat CId"
          cellRenderer={row => <Cell key={row+''}>{externalRecordsSet[row].variantCId}</Cell>}
        />
        <Column
          name="External Source Name"
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
          name="Distance To Feature"
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
      let externalRecords = variantExternalRecords.externalRecords;
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

