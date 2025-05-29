import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/term'
import { RegulatoryElementTerm } from '@/components/core/term'
import { VariantSet } from '@/components/core/set'
import { Table, Cell, Column} from '@/app/components/Table'
import { z } from 'zod'
import { getRegulatoryElementPosition } from '@/components/service/regulatoryElementInfo'
import { downloadBlob } from '@/utils/download'
import { resolveVariantCaID, variantIdResolveErrorMessage, linkedDataHubErroMessage } from './variantUtils'
import { getGitDataHubVariantInfo } from './variantInfoSources/gitDataHubVariantInfo'


export const GetRegulatoryElementsForThisVariant = MetaNode('GetRegulatoryElementForThisVariant')
  .meta({
    label: 'Identify regulatory elements associated with variant',
    description: 'Retrieve regulatory elements the overlaps with the given variant(s) from CFDE LDH.',
  })
  .inputs({ variant: VariantTerm })
  .output(RegulatoryElementTerm)
  .resolve(async (props) => {
    var varCaId = await resolveVariantCaID(props.inputs.variant);
    if(varCaId == null || varCaId == ''){
      throw new Error(variantIdResolveErrorMessage);
    }

    const response = await getGitDataHubVariantInfo(varCaId);
    if(response == null || response.data == null){
      throw new Error(linkedDataHubErroMessage);
    }

    if(response.data.ldFor.RegulatoryElement != null && response.data.ldFor.RegulatoryElement.length == 1){
      return response.data.ldFor.RegulatoryElement[0].entId;
    }
    return "N/A";
  }).story(props => ({
    abstract: `Regulatory element(s) that overlap with the given variant(s) were retrieved from the CFDE Llinked Data Hub API results\\ref{CFDE Linked Data Hub, https://ldh.genome.network/cfde/ldh/}`,
    introduction: `CFDE LDH is a graph-based network that facilitates access to excerpted regulatory information from external databases and studies including SCREEN, GTEx, and EN-TEx\\ref{doi:10.1038/s41586-020-2493-4}\\ref{doi:10.1126/science.aaz1776}\\ref{doi:10.1126/science.aar3146}`,
    methods: `Input variant(s) were queried through CFDE LDH API endpoints. Regulatory elements from SCREEN, eCLIP, and other databases that are linked to the variant(s) in LDH were retrieved from the JSON response.`,
    tableLegend: `A list of regulatory elements from CFDE LDH that overlap with the given variant(s)`,
  })).build()

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
      <>
        <p style={{fontSize: '14px'}}><b>Note:</b> In order to view all data, if avaliable, please expand the table rows!</p>
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
            name="Variant CAID"
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
      </> 
    )
  })
  .build()

  export const GetRegulatoryElementsForVariantSet = MetaNode('GetRegulatoryElementsForVariantSet')
  .meta({
    label: `Identify regulatory elements associated with variant(s)`,
    description: "Retrieve regulatory element(s) that overlap with each given variant(s) from CFDE LDH.",
    external: true,
  })
  .inputs({ variantset: VariantSet })
  .output(REforVariantSetInfo)
  .resolve(async (props) => {
    var variantSetInpt = props.inputs.variantset.set;
    
    let varAndRegulatoryElem = [];

    for(let indx in variantSetInpt){
      let varCaID = variantSetInpt[indx];
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

        const rePositionData = await getRegulatoryElementPosition(reEntId);
        let coordinates = null;
        if(rePositionData != null && rePositionData.coordinates != null){
          coordinates = rePositionData.coordinates;
        }else{
          continue;
        }
        tempObj.rePosition = coordinates.chromosome+":"+coordinates.start+"-"+coordinates.end;

        varAndRegulatoryElem.push(tempObj);
      }
    }

    if(varAndRegulatoryElem.length == 0){
      throw new Error(linkedDataHubErroMessage);
    }

    return varAndRegulatoryElem;
  }).story(props => ({
    abstract: `Regulatory element(s) that overlap with the given variant(s) were retrieved from the CFDE LDH API results\\ref{CFDE Linked Data Hub, https://ldh.genome.network/cfde/ldh/}`,
    introduction: `CFDE LDH is a graph-based network that facilitates access to excerpted regulatory information from external databases and studies including SCREEN, GTEx, and EN-TEx\\ref{doi:10.1038/s41586-020-2493-4}\\ref{doi:10.1126/science.aaz1776}\\ref{doi:10.1126/science.aar3146}`,
    methods: `Input variant(s) were queried through CFDE LDH API endpoints. Regulatory elements from SCREEN, eCLIP, and other databases that are linked to the variant(s) in LDH were retrieved from the JSON response.`,
    tableLegend: `A list of regulatory elements from CFDE LDH that overlap with the given variant(s)`,
  })).build()
