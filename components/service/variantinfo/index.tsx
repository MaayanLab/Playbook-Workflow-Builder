import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/input/term'
import { RegulatoryElementSet } from '@/components/core/input/set'
import { Table, Cell, Column, EditableCell } from '@/app/components/Table'
import { z } from 'zod'
//import { varinat_icon, variantinfo_icon } from '@/icons'

export const MyVariantInfoC =  z.object({
  data: z.object({
    entId: z.string(),
    entType: z.string(),
    ld: z.object({
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
            }).optional()
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

export type MyVariantInfo = z.infer<typeof MyVariantInfoC>

const MyAlleleRegistryExterSourcesListC = z.array(
  z.object({ 
    name: z.string(),
    sources: z.array(z.object({ 
      '@id':z.string(),
       id: z.string()
    }))
  })
);

export const VariantInfo = MetaNode('VariantInfo')
  .meta({
    label: 'Variant Iformation',
    description: 'A Variant resolved with MyVariantInfo',
    //icon: [varinat_icon]
  })
  .codec(MyVariantInfoC)
  .view(varinatinfo => (
    <div>
      <a target="_blank" href={`https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=${varinatinfo.data.entId}`}>{varinatinfo.data.entId}</a> (variant)
    </div>  
  ))
  .build()

async function myvariantinfo_query(variantId: string): Promise<MyVariantInfo> {
    const res = await fetch(`https://genboree.org/cfde-gene-dev/Variant/id/${encodeURIComponent(variantId)}`)
    return await res.json()
  }

export async function myAlleleInfo_query(variantId: string): Promise<object> {
    const res = await fetch(`https://reg.genome.network/allele/${encodeURIComponent(variantId)}`);
    return await res.json()
}

export const VarinatInfoFromVariantTerm = MetaNode('VarinatInfoFromVariantTerm')
  .meta({
    label: 'Resolve Variant Info from Term',
    description: 'Resolve variant info from variant term with MyVarinatInfo',
    //icon: [variantinfo_icon],
  })
  .inputs({ variant: VariantTerm })
  .output(VariantInfo)
  .resolve(async (props) => {
    return await myvariantinfo_query(props.inputs.variant);
  })
  .build()


export const GetRegulatoryElementsForVariantInfo = MetaNode('GetRegulatoryElementsForVariantInfo')
.meta({
  label: 'Resolve Reg. Elements from Var. Info',
  description: 'GetRegulatoryElementsForVariantInfo',
})
.inputs({ variantInfo: VariantInfo  })
.output(RegulatoryElementSet)
.resolve(async (props) => {
  return props.inputs.variantInfo.data.ldFor.RegulatoryElement.map(({ entId }) => entId)
})
.build()

export const xQTL_EvidenceDataTable = MetaNode('xQTL_EvidenceDataTable')
  .meta({
    label: 'xQTL_EvidenceDataTable',
    description: ''
  })
  .codec(MyVariantInfoC)
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
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].entContent.sourceDescription}</Cell>}
        />
        <Column
          name="Normalized Effect Size (nes)"
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].entContent.esQTL.nes}</Cell>}
        />
        <Column
          name="p-value"
          cellRenderer={row => <Cell key={row+''}>{xqtlEvidences[row].entContent.esQTL.sig}</Cell>}
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
    let xqtlEvidences = props.inputs.variantInfo.data.ld.xqtlEvidence;
      for(let e in xqtlEvidences){
        let xqtlE_entContent = xqtlEvidences[e].entContent;
        if(xqtlE_entContent.hasOwnProperty('eQTL')){
          xqtlE_entContent.esQTL = xqtlE_entContent.eQTL;
          delete xqtlE_entContent.eQTL;         
        }else if(xqtlE_entContent.hasOwnProperty('sQTL')){
          xqtlE_entContent.esQTL = xqtlE_entContent.sQTL;
          delete xqtlE_entContent.sQTL;
        }
      }
    return props.inputs.variantInfo;
  })
  .build()

  export const AlleleRegistryExternalSourcesTable = MetaNode('AlleleRegistryExternalSourcesTable')
  .meta({
    label: 'AlleleRegistryExternalSourcesTable',
    description: ''
  })
  .codec(MyAlleleRegistryExterSourcesListC)
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
              <Cell key={row+''}>
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

  export const GetAlleleRegistryExternalRecordsForVariantList = MetaNode('GetAlleleRegistryExternalRecordsForVariantList')
  .meta({
    label: 'Resolve Allele Registry External Records for Variant List',
    description: 'GetAlleleRegistryExternalRecordsForVariantList',
  })
  .inputs({ variantInfo: VariantInfo  })
  .output(AlleleRegistryExternalSourcesTable)
  .resolve(async (props) => {
    let alleleInfoExternalResources = [];
    let alleleInfo : any = await myAlleleInfo_query(props.inputs.variantInfo.data.entId);

    if(alleleInfo != null && alleleInfo['externalRecords'] != null){      
      let externalSources = alleleInfo['externalRecords'];
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

