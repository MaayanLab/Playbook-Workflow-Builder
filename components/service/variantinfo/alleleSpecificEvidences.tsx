import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/term'
import { VariantSet } from '@/components/core/set'
import { Table, Cell, Column} from '@/app/components/Table'
import { z } from 'zod'
import { downloadBlob } from '@/utils/download'
import { resolveVariantCaID, variantIdResolveErrorMessage, linkedDataHubErroMessage } from './variantUtils'
import { getVariantSetInfo } from './variantInfoSources/alleleRegistryVariantInfo'
import { getGitDataHubVariantInfo } from './variantInfoSources/gitDataHubVariantInfo'

export const AlleleSpecificEvidenceInfoC = z.array(
    z.object({
      ldhId: z.string(),
      ldhIri: z.string(),
      sourceDescription: z.string(),
      kbIri: z.string(),
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
        <>
          <p style={{fontSize: '14px'}}><b>Note:</b> In order to view all data, if avaliable, please expand the table rows!</p>
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
              name="Tissue Site or Cell Type"
              cellRenderer={row => <Cell key={row+''}>{alleleSpecificEvidence[row].sourceDescription.replace(/_/g, " ")}</Cell>}
            />
            <Column
              name="Evidence doc link"
              cellRenderer={row => <Cell key={row+''}><a target="_blank" href={`${alleleSpecificEvidence[row].ldhIri}`}>evidence link</a></Cell>}
            />
            <Column
              name="Source Iri"
              cellRenderer={row => <Cell key={row+''}><a target="_blank" href={`${alleleSpecificEvidence[row].kbIri}`}>Source Iri link</a></Cell>}
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
              name="Ref. allele level"
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
              name="Alt. allele level"
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
        </>  
      )
  })
  .build()

function getAlleleSpecificEvdncFromGitDataHub(alleleSpecificEvidencesList: any){
    var alleleSpecificEvidence: any = [];
    for(let a in alleleSpecificEvidencesList){
      let entContent = alleleSpecificEvidencesList[a].entContent;
      if(entContent == null){
        continue;
      }

      let kbIri = 'NA';
      if(entContent.kbIri != null){
        kbIri = entContent.kbIri;
      }

      let sourceDescription = 'NA';
      if(entContent.sourceDescription != null){
        sourceDescription = entContent.sourceDescription;
      }

      let specificities = entContent.AlleleSpecificity;
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
        "sourceDescription": sourceDescription,
        "kbIri": kbIri,
        "alleleSpecificityList": alleleSpecificityList
      }
      alleleSpecificEvidence.push(specificEvidence)
    }
    return alleleSpecificEvidence;
  }
  
  export const GetAlleleSpecificEvidencesForThisVariant = MetaNode('GetAlleleSpecificEvidencesForThisVariant')
    .meta({
      label: 'Identify associated allele specific epigenomic signature',
      description: 'Retrieve variant associated allele specific epigenomic signatures based on Roadmap and ENTEx data',
      external: true,
    })
    .inputs({ variant: VariantTerm })
    .output(AlleleSpecificEvidencesTable)
    .resolve(async (props) => {
      var varCaId = await resolveVariantCaID(props.inputs.variant);
      if(varCaId == null || varCaId == ''){
        throw new Error(variantIdResolveErrorMessage);
      }
  
      const response = await getGitDataHubVariantInfo(varCaId);
      if(response == null || response.data == null){
        throw new Error(linkedDataHubErroMessage);
      }
  
      let alleleSpecificEvidencesList = null;
      if(response.data.ld != null &&  response.data.ld.AlleleSpecificEvidence != null){
        alleleSpecificEvidencesList = response.data.ld.AlleleSpecificEvidence;
      }
      return getAlleleSpecificEvdncFromGitDataHub(alleleSpecificEvidencesList);
    }).story(props => ({
      abstract: `Variant/variant set associated allele specific epigenomic signatures were retrieved from CFDE LDH\\ref{Genomic Location Registry, https://reg.genome.network/reg/loc/} based on Roadmap and ENTEx data\\ref{doi:10.1038/nature14248}\\ref{doi:10.1126/science.aar3146}.`,
      introduction: `CFDE LDH is a graph-based network that facilitates access to excerpted regulatory information from external databases and studies including SCREEN, GTEx, and EN-TEx\\ref{doi:10.1038/s41586-020-2493-4}\\ref{doi:10.1126/science.aaz1776}\\ref{doi:10.1126/science.aar3146}`,
      methods: `Input variant(s) were queried through CFDE LDH API endpoints and their associated allele specificity evidence was retrieved in JSON response.`,
      tableLegend: `A table displaying 1) the type of allele specific epigenomic signature identified at the variant site, 2) the level of epigenomic modification at either allele, 3) its statistical significance, and 4) the original source of the evidence`,
    })).build()

    export const AlleleSpecificEvidenceForVariantSet = MetaNode('AlleleSpecificEvidenceForVariantSet')
    .meta({
      label: 'Allele Specific Evidence For Variant Set',
      description: ''
    })
    .codec(z.object({
        total: z.number(),
        found: z.number(),
        result: z.array(z.object({
          caid: z.string(),
          alleleSpecificEvidence: AlleleSpecificEvidenceInfoC
        }))
      })
    )
    .view( alleleEvidncForVarSetObj => {

      let alleleEvidncForVarSetArray = alleleEvidncForVarSetObj.result;

      return ( 
        <>  
          <p style={{fontSize: '14px'}}><b>Note:</b> In order to view all data, if avaliable, please expand the table rows!</p>
          <p style={{margin:'15px 0px 15px 0px'}}>Allele specificity evidence found for {alleleEvidncForVarSetObj.found} out of {alleleEvidncForVarSetObj.total} variants queried.</p>
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
              name="Variant CAid"
              cellRenderer={row => <Cell key={row+''}>{alleleEvidncForVarSetArray[row].caid}</Cell>}
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
              name="Evidence doc link"
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
              name="Source Iri"
              cellRenderer={row =>
              <Cell key={row+''}>{
                <table style={{borderCollapse: 'collapse', width:'100%'}}>
                {alleleEvidncForVarSetArray[row].alleleSpecificEvidence.map(sources =>
                    <tr><td><a target="_blank" href={`${sources.kbIri}`}>Source Iri link</a></td></tr>
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
              name="Ref. allele level"
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
              name="Alt. allele level"
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
        </>
      )
    })
    .build()

    export const GetVariantSetAlleleSpecificEvidence = MetaNode('GetVariantSetAlleleSpecificEvidence')
    .meta({
      label: `Identify Associated Allele Specific Epigenomic Signature`,
      description: "Retrieve variant/variant set associated allele specific epigenomic signatures based on Roadmap and ENTEx data.",
      external: true,
    })
    .inputs({ variantset: VariantSet})
    .output(AlleleSpecificEvidenceForVariantSet)
    .resolve(async (props) => {
      let variantSet = props.inputs.variantset.set;
      let variantSetInfo = await getVariantSetInfo(variantSet);
      if(variantSetInfo == null){
          throw new Error("No data available!");
      }
  
      let variantSetAlleleSpecificEvdnc = [];
      let cFound: number = 0;
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
        
        variantSetAlleleSpecificEvdnc.push(tempObj);
        cFound++
      }
  
      if(variantSetAlleleSpecificEvdnc.length == 0){
        throw new Error(linkedDataHubErroMessage);
      }
  
      let returnObj = {
        total: variantSetInfo.length,
        found: cFound,
        result: variantSetAlleleSpecificEvdnc
      }
      return returnObj;
    }).story(props => ({
      abstract: `Variant/variant set associated allele specific epigenomic signatures were retrieved from CFDE LDH\\ref{Genomic Location Registry, https://reg.genome.network/reg/loc/} based on Roadmap and ENTEx data\\ref{doi:10.1038/nature14248}, \\ref{doi:10.1126/science.aar3146}.`,
      introduction: `CFDE LDH is a graph-based network that facilitates access to excerpted regulatory information from external databases and studies including SCREEN, GTEx, and EN-TEx`,
      methods: `Input variant(s) were queried through CFDE LDH API endpoints and their associated tissue- and cell type-specific allele specificity evidence was retrieved from the JSON response.`,
      tableLegend: `A table displaying 1) the type of allele specific epigenomic signature identified at the variant site, 2) the level of epigenomic modification at either allele, 3) its statistical significance, and 4) the original source of the evidence`,
    })).build()