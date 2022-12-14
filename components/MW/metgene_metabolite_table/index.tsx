import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'


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
        <table>
          <tr>
            <th>{heading2}</th>
            <th>{heading3}</th>
            <th>{heading4}</th>
            <th>{heading5}</th>
          </tr>
          {data.map((val, key) => {
          var all_rxn_ids = val.KEGG_REACTION_ID
          var cpd_id = val.KEGG_COMPOUND_ID     
          var rxn_id_arr = all_rxn_ids.split(", ")
          
          function wrapURLRxn(rxn_id) {
            return "https://www.genome.jp/entry/rn:"+rxn_id  
          }
          var rxn_id_urls = rxn_id_arr.map(wrapURLRxn)
          
          var rxn_urls = []
          for (var i = 0; i < rxn_id_arr.length; i++){
            var rxn_id = rxn_id_arr[i]+" "
            var rxn_url = rxn_id_urls[i]
            
            rxn_urls.push(<a href = {rxn_url} target = "_blank">{rxn_id}</a>)
            
          }
          return (
            <tr key={key}>
              <td><a href = {`https://www.kegg.jp/entry/${val.KEGG_COMPOUND_ID}`} target = "__blank">{val.KEGG_COMPOUND_ID}</a></td>
              <td><a href = {`https://www.metabolomicsworkbench.org/databases/refmet/refmet_details.php?REFMET_NAME=${val.REFMET_NAME}`} target = "_blank">{val.REFMET_NAME}</a></td>
              <td>{rxn_urls}</td>
              <td><a href = {val.METSTAT_LINK} target = "__blank">METSTAT</a></td>
            </tr>
          )
        })}
        </table>
      </div>
    )
  })
  .build()
