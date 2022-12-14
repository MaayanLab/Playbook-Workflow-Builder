import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'

// A unique name for your data type is used here
export const MetGeneRxnTable = MetaNode.createData('MetGeneRxnTable')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGene Reaction Table',
    description: 'MetGene Reaction Table, rendered',
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using zod, a compile-time and runtime type-safe codec can be constructed
  .codec(z.any())
  
  // react component rendering your data goes here
  
  .view(data => {
    
    const heading1 = "Gene"
    const heading2 = "Kegg Rxn Id"
    const heading3 = "Kegg Rxn Name"
    
    return (
     
      <div>
        <h2>MetGENE Reactions</h2>
        <table>
          <tr>
            <th>{heading1}</th>
            <th>{heading2}</th>
            <th>{heading3}</th>
          </tr>
          {data.map((val, key) => {
          return (
            <tr key={key}>
              <td>{val.Gene}</td>
              <td>{val.KEGG_REACTION_ID}</td>
              <td>{val.KEGG_REACTION_NAME}</td>
            </tr>
          )
        })}
        </table>
      </div>
    )
  })
  .build()
