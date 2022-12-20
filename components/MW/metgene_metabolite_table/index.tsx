import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { SERVER_PROPS_ID } from 'next/dist/shared/lib/constants'


// A unique name for your data type is used here
export const MetgeneMetaboliteTable = MetaNode.createData('MetgeneMetaboliteTable')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'Metgene metabolite table',
    description: 'Metgene metabolite table, rendered',
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using zod, a compile-time and runtime type-safe codec can be constructed
  .codec(z.any())
  // react component rendering your data goes here
  .view(data => {
    const heading1 = "Gene"
    const heading2 = "KEGG_COMPOUND_ID"
    const heading3 = "REFMET_NAME"
    const heading4 = "KEGG_REACTION_ID"
    const heading5 = "METSTAT_LINK"
    
    return (
      <div>
        <h2>MetGENE metabolites</h2>
        {data.map((geneval:any, index:number) => (
          <div key={index}>
            
            <table>
              <thead>
              <tr>
                <th>{heading1}</th>
                <th>{heading2}</th>
                <th>{heading3}</th>
                <th>{heading4}</th>
                <th>{heading5}</th>
              </tr>
              </thead>
              <tbody>
              
                {geneval.map((val:any, index:number) => {
                  const all_rxn_ids = val.KEGG_REACTION_ID
                  const cpd_id = val.KEGG_COMPOUND_ID 
                  const rxn_id_arr = all_rxn_ids.split(", ")
                            
                  
                    return (
                      <tr key={index}>
                        <td>{val.Gene}</td>
                        <td><a href = {`https://www.kegg.jp/entry/${val.KEGG_COMPOUND_ID}`} target = "__blank">{val.KEGG_COMPOUND_ID}</a></td>
                        <td><a href = {`https://www.metabolomicsworkbench.org/databases/refmet/refmet_details.php?REFMET_NAME=${val.REFMET_NAME}`} target = "_blank">{val.REFMET_NAME}</a></td>
                        <td>
                          {rxn_id_arr.map((rxn_id:any, i:number) => {
                            let sep_char = ", ";
                            if (i == (rxn_id_arr.length-1)) {
                              sep_char = "";
                            }
                            return (
                            <a key={i} href = {`https://www.genome.jp/entry/rn:${rxn_id}`} target = "_blank">{rxn_id}{sep_char}</a>
                          
                          )})}
                        </td>
                        <td><a href = {val.METSTAT_LINK} target = "__blank">METSTAT</a></td>
                      </tr>
                    )
                })}
              
              </tbody>
            </table>
          </div>
        ))}

      </div>

    )
  })
  .build()
