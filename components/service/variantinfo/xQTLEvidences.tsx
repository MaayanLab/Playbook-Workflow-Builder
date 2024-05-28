import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/term'
import { VariantSet } from '@/components/core/set'
import { Table, Cell, Column} from '@/app/components/Table'
import { z } from 'zod'
import { downloadBlob } from '@/utils/download'
import { resolveVariantCaID, variantIdResolveErrorMessage, gitDataHubErroMessage } from './variantUtils'
import { getVariantSetInfo } from './variantInfoSources/alleleRegistryVariantInfo'
import { getGitDataHubVariantInfo, GitHubVariantInfoC } from './variantInfoSources/gitDataHubVariantInfo'
import { xQTL_EvidenceArray } from './variantInfoSources/gitDataHubVariantInfo'

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
    var varCaId = await resolveVariantCaID(props.inputs.variant);
    if(varCaId == null || varCaId == ''){
      throw new Error(variantIdResolveErrorMessage);
    }

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
          name="Variant CAID"
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

  export const GetVariantSetXQTLEvidence = MetaNode('GetVariantSetXQTLEvidence')
  .meta({
    label: `Get Variant Set xQTL Evidence`,
    description: "Get xQTL Evidence form Variant Set."
  })
  .inputs({ variantset: VariantSet })
  .output(xQTL_EvidenceFroVariantSet)
  .resolve(async (props) => {
    let variantSet = props.inputs.variantset.set;
    let variantSetInfo = await getVariantSetInfo(variantSet);
    if(variantSetInfo == null){
        throw new Error("No data available!");
    }

    let variantSetXQTLEvidnc = [];
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

      variantSetXQTLEvidnc.push(tempObj);
    }

    if(variantSetXQTLEvidnc.length == 0){
      throw new Error(gitDataHubErroMessage);
    }

    return variantSetXQTLEvidnc;
  }).story(props =>
    "Get xQTL Evidence form Variant Set."
  ).build()
