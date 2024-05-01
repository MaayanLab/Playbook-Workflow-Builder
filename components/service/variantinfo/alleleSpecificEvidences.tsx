import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/term'
import { Table, Cell, Column} from '@/app/components/Table'
import { z } from 'zod'
import { downloadBlob } from '@/utils/download'
import { resolveVarinatCaID, variantIdResolveErrorMessage, gitDataHubErroMessage } from './variantUtils'
import { VariantSetInfo } from './variantInfoSources/alleleRegistryVariantInfo'
import { getGitDataHubVariantInfo } from './variantInfoSources/gitDataHubVariantInfo'

export const AlleleSpecificEvidenceInfoC = z.array(
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