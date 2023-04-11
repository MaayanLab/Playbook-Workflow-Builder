import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/input/term'
import { RegulatoryElementTerm } from '@/components/core/input/term'
import { Table, Cell, Column, EditableCell } from '@/app/components/Table'
import { z } from 'zod'
//import { varinat_icon, variantinfo_icon } from '@/icons'

const AlleleRegistryVariantInfoC = z.object({
  '@id': z.string(),
  entId: z.string(),
  externalRecords: z.any({}).optional()
})
export type AlleleRegistryVariantInfo = z.infer<typeof AlleleRegistryVariantInfoC>

const AlleleRegistryExternalSourcesInfoC = z.array(
  z.object({ 
    name: z.string(),
    sources: z.array(z.object({ 
      '@id':z.string(),
      id: z.string()
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
            AlleleSpecificity:z.any({}).optional()
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

export const VariantInfo = MetaNode('VariantInfo')
  .meta({
    label: 'Variant Iformation',
    description: 'A Variant resolved with MyVariantInfo',
    //icon: [varinat_icon]
  })
  .codec(AlleleRegistryVariantInfoC)
  .view(varinatinfo => (
    <div>
      <a target="_blank" href={`https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=${varinatinfo.entId}`}>{varinatinfo.entId}</a> (variant)
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

export const VarinatInfoFromVariantTerm = MetaNode('VarinatInfoFromVariantTerm')
  .meta({
    label: 'Resolve Variant Info from Term',
    description: 'Resolve variant info (Allele registry API) from variant term with MyVarinatInfo',
    //icon: [variantinfo_icon],
  })
  .inputs({ variant: VariantTerm })
  .output(VariantInfo)
  .resolve(async (props) => {
    const response = await getAlleleRegistryVariantInfo(props.inputs.variant);
    response.entId = props.inputs.variant;
    return response;
  })
  .build()


export const GetRegulatoryElementsForThisVariant = MetaNode('GetRegulatoryElementsForThisVariant')
.meta({
  label: 'Resolve Reg. Elements from Var. Info',
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
.codec(GitHubVariantInfoC)
.view(varinatinfo => {
    let alleleSpecificEvidences = varinatinfo.data.ld.AlleleSpecificEvidence;

      return (
          <Table
            height={500}
            cellRendererDependencies={[alleleSpecificEvidences]}
            numRows={alleleSpecificEvidences.length}
            enableGhostCells
            enableFocusedCell
          >
            <Column
              name="ldhId"
              cellRenderer={row => <Cell key={row+''}>{alleleSpecificEvidences[row].ldhId}</Cell>}
            />
            <Column
              name="ldhIri"
              cellRenderer={row => <Cell key={row+''}><a target="_blank" href={`${alleleSpecificEvidences[row].ldhIri}`}>evidence link</a></Cell>}
            />
            <Column
              name="Source description"
              cellRenderer={row => <Cell key={row+''}>{alleleSpecificEvidences[row].entContent.sourceDescription.replace(/_/g, " ")}</Cell>}
            />
          </Table>
      )
}).build()
  /*
  <Column
    name="Reference Allele Quant"
    cellRenderer={row => <Cell key={row+''}>{alleleSpecificEvidences[row].entContent.AlleleSpecificity.refAlleleQuant}</Cell>}
  />
  <Column
    name="Alternate Allele Quant"
    cellRenderer={row => <Cell key={row+''}>{alleleSpecificEvidences[row].entContent.AlleleSpecificity.altAlleleQuant}</Cell>}
  />*/

export const GetAlleleSpecificEvidencesForThisVariant = MetaNode('GetAlleleSpecificEvidencesForThisVariant')
.meta({
  label: 'Resolve Allele Specific Evidences from Var. Info',
  description: 'GetAlleleSpecificEvidencesForThisVariant',
})
.inputs({ variantInfo: VariantInfo  })
.output(AlleleSpecificEvidencesTable)
.resolve(async (props) => {
  return await getGitDataHubVariantInfo(props.inputs.variantInfo.entId);
})
.build()

export const xQTL_EvidenceDataTable = MetaNode('xQTL_EvidenceDataTable')
  .meta({
    label: 'xQTL_EvidenceDataTable',
    description: ''
  })
  .codec(GitHubVariantInfoC)
  .view(varinatinfo => {
    let xqtlEvidences = varinatinfo.data.ld.xqtlEvidence;
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
          name="Tissue site"
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].entContent.sourceDescription.replace(/_/g, " ")}</Cell>}
        />
        <Column
          name="Normalized Effect Size (nes)"
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].entContent.esQTL.nes}</Cell>}
        />
        <Column
          name="p-value"
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].entContent.esQTL.sig}</Cell>}
        />
        <Column
          name="Type"
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].entContent.type}</Cell>}
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
              name="Source Id and Link"            
              cellRenderer={row => 
              <Cell  key={row+''}>
                  {AlleleRegistryExternalSourcesList[row].sources.map(sources =>
                    <table style={{borderCollapse: 'collapse', width:'100%'}}>
                      <tr>
                        <td style={{width:'50%'}}>{ sources.id }</td>
                        <td style={{width:'50%'}}><a target="_blank" href={`${sources['@id']}`}>Resource link</a></td>             
                      </tr>
                    </table>
                  )}
              </Cell>}
            />
          </Table>
        )
}).build()

  export const GetAlleleRegistryExternalRecordsForVariant = MetaNode('GetAlleleRegistryExternalRecordsForVariant')
  .meta({
    label: 'Resolve Allele Registry External Records for Variant',
    description: 'GetAlleleRegistryExternalRecordsForVariant',
  })
  .inputs({ variantInfo: VariantInfo  })
  .output(AlleleRegistryExternalRecordsTable)
  .resolve(async (props) => {
    let alleleInfoExternalResources = [];
    let variantInfoObj: any = props.inputs.variantInfo;

    if(variantInfoObj['externalRecords'] != null){      
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
    }   
    return alleleInfoExternalResources;
  })
  .build()

