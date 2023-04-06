import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/input/term'
import { Table, Cell, Column, EditableCell } from '@/app/components/Table'
import { z } from 'zod'
//import { varinat_icon, variantinfo_icon } from '@/icons'

export const MyVariantInfoC =  z.object({
    data: z.object({
      entId: z.string(),
      entType: z.string(),
      ld: z.object({
        xqtlEvidence: z.array(
            z.object({ entId: z.string() })
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

