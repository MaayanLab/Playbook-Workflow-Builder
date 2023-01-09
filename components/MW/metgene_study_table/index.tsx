import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'

const MetGeneStudyObjC = z.object({
  KEGG_COMPOUND_ID: z.string(),
  REFMET_NAME: z.string(),
  STUDY_ID: z.string(),
});

export type  MetGeneStudyObj = z.infer<typeof  MetGeneStudyObjC>

// important ref: https://rsinohara.github.io/json-to-zod-react/
// array of objects: [{}]
export const MetGeneStudyObjArrayC = z.array(
  MetGeneStudyObjC
)

export type MetGeneStudyObjArray = z.infer<typeof MetGeneStudyObjArrayC>

// A unique name for your data type is used here
export const MetGeneStudyTable = MetaNode.createData('MetGeneStudyTable')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGene Studies Table',
    description: 'Studies table corresponding to gene',
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using zod, a compile-time and runtime type-safe codec can be constructed
  .codec(MetGeneStudyObjArrayC)
  // react component rendering your data goes here
  .view(data => {
    const heading1 = "Kegg Compound Id"
    const heading2 = "Refmet Name"
    const heading3 = "Study Ids"
    
    const studyID_url = "https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID="
    return (
      <div>
        <h2>MetGENE Studies</h2>
        <table>
          <tr>
            <th>{heading1}</th>
            <th>{heading2}</th>
            <th>{heading3}</th>
          </tr>
          {data.map((val : MetGeneStudyObj, key : number) => {
          var all_study_ids = val.STUDY_ID
                
          var study_id_arr = all_study_ids.split(", ")
          
          
          return (
            <tr key={key}>
              <td>{val.KEGG_COMPOUND_ID}</td>
              <td><a href = {`https://www.metabolomicsworkbench.org/databases/refmet/refmet_details.php?REFMET_NAME=${val.REFMET_NAME}`} target = "_blank">{val.REFMET_NAME}</a></td>
              <td>
                {study_id_arr.map((study_id: string, i: number) => 
                  <>
                  {i > 0 ? <span>, </span> : null}
                  <a key={i} href = {`https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=${study_id}`} target = "_blank">{study_id}</a>
                  </>
                )}
              </td>
            </tr>
          )
        })}
        </table>
      </div>
    )
  })
  .build()
