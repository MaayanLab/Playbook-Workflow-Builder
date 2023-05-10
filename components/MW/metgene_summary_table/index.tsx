import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'

const MetGeneSummaryObjC = z.object({
  Pathways: z.number(),
  Reactions: z.number(),
  Metabolites: z.number(),
  Studies: z.number(),
  Genes: z.string(),
});

export type  MetGeneSummaryObj = z.infer<typeof  MetGeneSummaryObjC>

// important ref: https://rsinohara.github.io/json-to-zod-react/
// array of objects: [{}]
export const MetGeneSummaryObjArrayC = z.array(
  MetGeneSummaryObjC
)

export type MetGeneSummaryObjArray = z.infer<typeof MetGeneSummaryObjArrayC>
// A unique name for your data type is used here
export const MetGeneSummaryTable = MetaNode('MetGeneSummaryTable')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'Metgene Summary Table',
    description: 'Metgene Summary Table comprises of the number of pathways, reactions, metabolites and studies corresponding to the given gene',
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using zod, a compile-time and runtime type-safe codec can be constructed
  .codec(MetGeneSummaryObjArrayC)

  // react component rendering your data goes here
  .view(data => {
    const heading1 = "Gene"
    const heading2 = "Pathways"
    const heading3 = "Reactions"
    const heading4 = "Metabolites"
    const heading5 = "Studies"
    return (
      <div className="prose max-w-none">
        <h2>MetGENE Summary</h2>
        <table>
          <tr>
            <th>{heading1}</th>
            <th>{heading2}</th>
            <th>{heading3}</th>
            <th>{heading4}</th>
            <th>{heading5}</th>
          </tr>
          {data.map((val : MetGeneSummaryObj, key : number) => {



          return (
            <tr key={key}>
              <td><a href = {`https://bdcw.org/MetGENE/geneInfo.php?species=hsa&GeneIDType=SYMBOL&anatomy=NA&disease=NA&phenotype=NA&GeneInfoStr=${val.Genes}`} target = "_blank">{val.Genes}</a></td>
              <td><a href = {`https://bdcw.org/MetGENE/pathways.php?species=hsa&GeneIDType=SYMBOL&anatomy=NA&disease=NA&phenotype=NA&GeneInfoStr=${val.Genes}`} target = "_blank">{val.Pathways}</a></td>
              <td><a href = {`https://bdcw.org/MetGENE/reactions.php?species=hsa&GeneIDType=SYMBOL&anatomy=NA&disease=NA&phenotype=NA&GeneInfoStr=${val.Genes}`} target = "_blank">{val.Reactions}</a></td>
              <td><a href = {`https://bdcw.org/MetGENE/metabolites.php?species=hsa&GeneIDType=SYMBOL&anatomy=NA&disease=NA&phenotype=NA&GeneInfoStr=${val.Genes}`} target = "_blank">{val.Metabolites}</a></td>
              <td><a href = {`https://bdcw.org/MetGENE/studies.php?species=hsa&GeneIDType=SYMBOL&anatomy=NA&disease=NA&phenotype=NA&GeneInfoStr=${val.Genes}`} target = "_blank">{val.Studies}</a></td>
            </tr>
          )
        })}
        </table>
      </div>
    )
  })
  .build()
