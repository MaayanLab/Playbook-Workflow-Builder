import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/term'
import { Table, Cell, Column} from '@/app/components/Table'
import { z } from 'zod'
import { downloadBlob } from '@/utils/download'
import { resolveVarinatCaID, variantIdResolveErrorMessage, alleleRegRespErrorMessage } from './variantUtils'
import { getAlleleRegistryVariantInfo,  AlleleRegistryVariantInfo, VariantSetInfo } from './variantInfoSources/alleleRegistryVariantInfo'

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
          name="Variant CAID"
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
