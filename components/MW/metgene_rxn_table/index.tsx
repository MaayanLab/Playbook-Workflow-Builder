import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'

// object: {}
const MetGeneRxnObjC = z.object({
  Gene: z.string(),
  KEGG_REACTION_ID: z.string(),
  KEGG_REACTION_NAME: z.string(),
});

export type  MetGeneRxnObj = z.infer<typeof  MetGeneRxnObjC>

// important ref: https://rsinohara.github.io/json-to-zod-react/
// array of objects: [{}]
export const MetGeneRxnObjArrayC = z.array(
  MetGeneRxnObjC
)

export type MetGeneRxnObjArray = z.infer<typeof MetGeneRxnObjArrayC>

// array of array of objects: [[{}]]
export const MetGeneRxnObjArray2C = z.array(
  MetGeneRxnObjArrayC
)

export type MetGeneRxnObjArray2 = z.infer<typeof MetGeneRxnObjArray2C>

// A unique name for your data type is used here
export const MetGeneRxnTable = MetaNode('MetGeneRxnTable')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGENE Reaction Table',
    description: 'MetGENE Reaction Table, rendered',
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using zod, a compile-time and runtime type-safe codec can be constructed
  .codec(MetGeneRxnObjArray2C)
  
  // react component rendering your data goes here
  
  .view(data => {
    
    const heading1 = "Gene"
    const heading2 = "Kegg Rxn Id"
    const heading3 = "Kegg Rxn Name"
    
    return (
      
      <div>
        <h2>MetGENE Reactions</h2>
        {data.map((arrayVal:MetGeneRxnObjArray, index:number) => (
          <div key={index}>
          <table>
            
            <thead>
            <tr>
              <th>{heading1}</th>
              <th>{heading2}</th>
              <th>{heading3}</th>
            </tr>
            </thead>
             
            <tbody>    
            {arrayVal.map((val:MetGeneRxnObj, i:number) => {
            
              return (
                <tr key={i}>
                  <td>{val.Gene}</td>
                  <td><a href = {`https://www.kegg.jp/entry/rn:${val.KEGG_REACTION_ID}`} target = "__blank">{val.KEGG_REACTION_ID}</a></td>
                  <td>{val.KEGG_REACTION_NAME}</td>
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
