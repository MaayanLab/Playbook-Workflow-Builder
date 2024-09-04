import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/term'
import { VariantSet } from '@/components/core/set'
import { Table, Cell, Column} from '@/app/components/Table'
import { z } from 'zod'
import { downloadBlob } from '@/utils/download'
import { resolveVariantCaID, variantIdResolveErrorMessage, alleleRegRespErrorMessage } from './variantUtils'
import { getAlleleRegistryVariantInfo,  AlleleRegistryVariantInfo, getVariantSetInfo } from './variantInfoSources/alleleRegistryVariantInfo'

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

export const AlleleRegistryExternalRecordsTable = MetaNode('AlleleRegistryExternalRecordsTable')
  .meta({
    label: 'Allele Registry External Records Table',
    description: ''
  })
  .codec(AlleleRegistryExternalSourcesInfoC)
  .view(AlleleRegistryExternalSourcesList => {
    let sourcesList = AlleleRegistryExternalSourcesList ?? [];

    return (
      <>
        <p style={{fontSize: '14px'}}><b>Note:</b> In order to view all data, if avaliable, please expand the table rows!</p>
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
            name="Database Name"
            cellRenderer={row => <Cell key={row+''}>{sourcesList[row].name}</Cell>}
          />
          <Column
            name="Variant id"
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
            name="Source link"
            cellRenderer={row =>
            <Cell  key={row+''}>
                <table style={{borderCollapse: 'collapse', width:'100%'}}>
                      {sourcesList[row].sources.map(sources =>
                          <tr><td><a target="_blank" href={`${sources['@id']}`}>Source link</a></td></tr>
                      )}
                </table>
            </Cell>}
          />
        </Table>
      </>
    )
  })
  .build()

  function processExternalRecords(variantInfoObj: AlleleRegistryVariantInfo, caid: string){
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

    let arExternalR = {
      name: 'Allele Registry',
      sources: [{
        '@id': variantInfoObj['@id'], 
        id: caid+"("+variantInfoObj.communityStandardTitle[0]+")"
      }]
    }
    alleleInfoExternalResources.push(arExternalR);

    return alleleInfoExternalResources;
  }

  export const GetAlleleRegistryExternalRecordsForVariant = MetaNode('GetAlleleRegistryExternalRecordsForVariant')
  .meta({
    label: 'Retrieve alternative identifiers for variant',
    description: 'Retrieve MyVariant.info, dbSNP, gnomAD, and other common identifiers for given variant from ClinGen Allele Registry.',
  })
  .inputs({ variant: VariantTerm  })
  .output(AlleleRegistryExternalRecordsTable)
  .resolve(async (props) => {
    var varCaId = await resolveVariantCaID(props.inputs.variant);
    if(varCaId == null || varCaId == ''){
      throw new Error(variantIdResolveErrorMessage);
    }

    let variantInfoObj: any = await getAlleleRegistryVariantInfo(varCaId);
    if(variantInfoObj == null){
      throw new Error(alleleRegRespErrorMessage);
    }

    let reponse = null;
    if(variantInfoObj['externalRecords'] != null){
      reponse = processExternalRecords(variantInfoObj, varCaId);
    }
    return reponse;
  }).story(props => ({
    abstract: `MyVariant.info, dbSNP, gnomAD, and other common identifiers for the given variant(s) were retrieved from ClinGen Allele Registry\\ref{doi:10.1002/humu.23637}.`,
    introduction: `ClinGen Allele Registry is a variant naming and registration service that provide stable and unique identifiers for any possible alleles in the human genome\\ref{doi:10.1002/humu.23637}.`,
    methods: `Input variant id(s) were queried through ClinGen Allele Registry API endpoints corresponding to the id type and alternative ids and aliases were retreived from the JSON response.`,
    legend: `A table displaying all the commonly used identifiers for the given variant(s)/allele(s)`,
  })).build()

  export const VariantSetExternalRecordsInfo = MetaNode('VariantSetExternalRecordsInfo')
  .meta({
    label: 'Variant Set External Records Info',
    description: ''
  })
  .codec(AlleleRegistryExternalSourcesSetInfoC)
  .view( externalRecordsSet => {
    return ( 
      <>
        <p style={{fontSize: '14px'}}><b>Note:</b> In order to view all data, if avaliable, please expand the table rows!</p>
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
            name="Variant CAid"
            cellRenderer={row => <Cell key={row+''}>{externalRecordsSet[row].variantCaId}</Cell>}
          />
          <Column
            name="Database name"
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
            name="Variant id"
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
            name="Source link "
            cellRenderer={row =>
              <Cell key={row+''}>
                <table style={{borderCollapse: 'collapse', width:'100%'}}>
                  {externalRecordsSet[row].externalRecords?.map(externalRecord =>
                      <tr><td><a target="_blank" href={externalRecord.sources[0]['@id']}>Source link</a></td></tr>
                  )}
                </table>
              </Cell>}
            />
        </Table>
      </>
    )  
  })
  .build()

  export function getExternalRecordsFromAlleleRegistry(variantSetInfo: any){
    let externalRecordsSet = [];
    for(let indx in variantSetInfo){
      let variantInfoObj = variantSetInfo[indx];
      if(variantInfoObj['externalRecords'] != null){
        let exteralRecordsData = processExternalRecords(variantInfoObj, variantInfoObj.entId);
        let variantExternalRecords = {
          "variantCaId" : variantInfoObj.entId,
          "externalRecords" : exteralRecordsData
        };

        externalRecordsSet.push(variantExternalRecords);
      }
    }
    return externalRecordsSet;
  }


  export const GetVariantSetExternalRecords = MetaNode('GetVariantSetExternalRecords')
  .meta({
    label: `Retrieve alternative identifiers for variant(s)`,
    description: "Retrieve MyVariant.info, dbSNP, gnomAD, and other common identifiers for given variant(s) from ClinGen Allele Registry."
  })
  .inputs({ variantset: VariantSet })
  .output(VariantSetExternalRecordsInfo)
  .resolve(async (props) => {
    let variantSet = props.inputs.variantset.set;
    let variantSetInfo = await getVariantSetInfo(variantSet);
    if(variantSetInfo == null){
        throw new Error("No data available for the inputed variant set! "+variantIdResolveErrorMessage);
    }

    return getExternalRecordsFromAlleleRegistry(variantSetInfo);
  }).story(props => ({
    abstract: `MyVariant.info, dbSNP, gnomAD, and other common identifiers for the given variant(s) were retrieved from ClinGen Allele Registry\\ref{doi:10.1002/humu.23637}.`,
    introduction: `ClinGen Allele Registry is a variant naming and registration service that provide stable and unique identifiers for any possible alleles in the human genome\\ref{doi:10.1002/humu.23637}.`,
    methods: `Input variant id(s) were queried through ClinGen Allele Registry API endpoints corresponding to the id type and alternative ids and aliases were retreived from the JSON response.`,
    legend: `A table displaying all the commonly used identifiers for the given variant(s)/allele(s)`,
  })).build()
