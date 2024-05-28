import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/term'
import { VariantSet } from '@/components/core/set'
import { Table, Cell, Column} from '@/app/components/Table'
import { z } from 'zod'
import { downloadBlob } from '@/utils/download'
import { resolveVariantCaID, variantIdResolveErrorMessage, getMyVarintInfoLink, alleleRegRespErrorMessage } from './variantUtils'
import { getAlleleRegistryVariantInfo, getVariantSetInfo } from './variantInfoSources/alleleRegistryVariantInfo'
import { AlleleRegistryExternalRecordsTable, VariantSetExternalRecordsInfo, getExternalRecordsFromAlleleRegistry } from './externalRecords'
import { getVariantInfoFromMyVariantInfo } from './variantInfoSources/myVariantInfo'

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
    if(Array.isArray(associatedGeens)){
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
    }else{
      associatedGeensMap[associatedGeens.gene_id+""] = associatedGeens;
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

  export const GetGeneForVariantFromMyVariantInfo = MetaNode('GetGeneForVariantFromMyVariantInfo')
  .meta({
    label: 'Identify Variant To Gene Association HG38',
    description: 'Get Associated Gene info for a given Variant. Use any of the most common types of Variant identifiers as an input value!'
  })
  .inputs({ variant: VariantTerm })
  .output(GeneAssociations_HG38)
  .resolve(async (props) => {
    var varCaId = await resolveVariantCaID(props.inputs.variant);
    if(varCaId == null || varCaId == ''){
      throw new Error(variantIdResolveErrorMessage);
    }

    const alleleRegResponse = await getAlleleRegistryVariantInfo(varCaId);
    if(alleleRegResponse == null){
      throw new Error(alleleRegRespErrorMessage);
    }

    let response = null;
    let myVariantInfoURL = getMyVarintInfoLink(alleleRegResponse);
    if(myVariantInfoURL != null){
      response = await getVariantInfoFromMyVariantInfo(myVariantInfoURL);
    }

    return processHG38ExternalRecordsResponse(response);
  })
  .story(props => `The closest gene to the variant was extract from the MyVariant.info API results [\\ref{doi:10.1093/bioinformatics/btac017}].`)
  .build()

  export const GetVariantToGeneAssociation_HG38 = MetaNode('GetVariantToGeneAssociation_HG38')
  .meta({
    label: `Identify Variant To Gene Association HG38`,
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
      response = await getVariantInfoFromMyVariantInfo(hg38ExtSourceLink);
    }

    return processHG38ExternalRecordsResponse(response);
  }).story(props =>
    `Get Associated Gene info for a given Variant.`
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
          name="Variant CAID"
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

  async function getGeneAssociationsHG38FromExternalRecords(variantExternalRecordsSetInfo: any){
    let variantGeneAssociationsSet = [];
    for(let vERIdx in variantExternalRecordsSetInfo){
      let variantExternalRecords = variantExternalRecordsSetInfo[vERIdx];
      let externalRecords: any = variantExternalRecords.externalRecords;
      for(let erIdx in externalRecords){
        if(externalRecords[erIdx] != null && externalRecords[erIdx].name == "MyVariantInfo_hg38"){
          let hg38ExtSourceLink = processHG38externalRecords(externalRecords[erIdx]);
          let response = await getVariantInfoFromMyVariantInfo(hg38ExtSourceLink);
          let associatedGeensList = processHG38ExternalRecordsResponse(response);

          let variantGeneAssociation = {
            "variantCaId": variantExternalRecords.variantCaId,
            "geneAssociationsHg38": associatedGeensList
          }
          variantGeneAssociationsSet.push(variantGeneAssociation);
        }
      }   
    }
    return variantGeneAssociationsSet;
  }

  export const GetVariantSetToGeneAssociation_HG38 = MetaNode('GetVariantSetToGeneAssociation_HG38')
  .meta({
    label: `Identify Variant (Set) And Gene Associations HG38`,
    description: "Get Associated Gene info for a given Variant Set."
  })
  .inputs({ variantset: VariantSet })
  .output(GeneAssociationsSet_HG38)
  .resolve(async (props) => {
    let variantSet = props.inputs.variantset.set;
    let variantSetInfo = await getVariantSetInfo(variantSet);
    if(variantSetInfo == null){
        throw new Error("No data available!");
    }

    let variantExternalRecordsSetInfo = getExternalRecordsFromAlleleRegistry(variantSetInfo);
    if(variantExternalRecordsSetInfo == null){
      throw new Error("No data avaliable!");
    }

    return await getGeneAssociationsHG38FromExternalRecords(variantExternalRecordsSetInfo);
  }).story(props =>
    `Get Associated Gene info for a given set of Variant External records.`
  ).build()


  export const GetVariantSetExternalRecToGeneAssociation_HG38 = MetaNode('GetVariantSetExternalRecToGeneAssociation_HG38')
  .meta({
    label: `Identify Variant (Set) And Gene Associations HG38`,
    description: "Get Associated Gene info for a given set of Variant External records."
  })
  .inputs({ variantSetExternalRecordsInfo: VariantSetExternalRecordsInfo })
  .output(GeneAssociationsSet_HG38)
  .resolve(async (props) => {
    let variantExternalRecordsSetInfo = props.inputs.variantSetExternalRecordsInfo;

    return await getGeneAssociationsHG38FromExternalRecords(variantExternalRecordsSetInfo);
  }).story(props =>
    `Get Associated Gene info for a given set of Variant External records.`
  ).build()
