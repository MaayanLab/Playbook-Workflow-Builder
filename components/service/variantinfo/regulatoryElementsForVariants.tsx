import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/term'
import { RegulatoryElementTerm } from '@/components/core/term'
import { VariantSet } from '@/components/core/set'
import { Table, Cell, Column} from '@/app/components/Table'
import { z } from 'zod'
import { getRegElemPositionData } from '@/components/service/regulatoryElementInfo'
import { downloadBlob } from '@/utils/download'
import { resolveVariantCaID, variantIdResolveErrorMessage, gitDataHubErroMessage } from './variantUtils'
import { getGitDataHubVariantInfo } from './variantInfoSources/gitDataHubVariantInfo'


export const GetRegulatoryElementsForThisVariant = MetaNode('GetRegulatoryElementForThisVariant')
  .meta({
    label: 'Identify Variant And Regulatory Element Association',
    description: 'Get Regulatory Elements For This Variant ID',
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
      throw new Error(gitDataHubErroMessage);
    }

    if(response.data.ldFor.RegulatoryElement != null && response.data.ldFor.RegulatoryElement.length == 1){
      return response.data.ldFor.RegulatoryElement[0].entId;
    }
    return "N/A";
  })
  .story(props => ({}))
  .build()

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
    label: `Identify regulatory elements associated with variant`,
    description: "Description change: Identify regulatory elements in the region of the variant."
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
    `Description change: Identify regulatory elements in the region of the variant.`
  ).build()
