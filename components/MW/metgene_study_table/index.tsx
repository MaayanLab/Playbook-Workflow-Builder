import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { additional_info_icon, study_icon } from '@/icons';

const MetGeneStudyObjC = z.object({
  KEGG_COMPOUND_ID: z.string().url(),
  REFMET_NAME: z.string(),
  STUDY_ID: z.string(),
});

export type MetGeneStudyObj = z.infer<typeof  MetGeneStudyObjC>

// important ref: https://rsinohara.github.io/json-to-zod-react/
// array of objects: [{}]
export const MetGeneStudyObjArrayC = z.array(
  MetGeneStudyObjC
)

export type MetGeneStudyObjArray = z.infer<typeof MetGeneStudyObjArrayC>

// A unique name for your data type is used here
export const MetGeneStudyTable = MetaNode('MetGeneStudyTable')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGENE Studies Table',
    description: 'Studies table corresponding to gene',
    icon: [study_icon, additional_info_icon],
    external: true,
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

    const contents = data.filter(arrval => arrval.STUDY_ID.split(', ').length > 0)
    if (contents.length === 0) {
      return (
        <div className="prose max-w-none">
          <h2 className="m-0">MetGENE Studies</h2>
          <span className="text-red-500">No studies found.</span>
        </div>
      )
    }

    return (
      <div className="prose max-w-none">
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
          const m = /href\s*=\s*["'](.+)["'].*?>?([^><]+)/g.exec(val.KEGG_COMPOUND_ID)
          const { url, compound } = m === null ?
            { url: `https://www.kegg.jp/entry/${val.KEGG_COMPOUND_ID}`, compound: val.KEGG_COMPOUND_ID }
            : { url: m[1], compound: m[2] }

          return (
            <tr key={key}>
              <>
              <td><a href = {url} target = "_blank">{compound}</a></td>
              <td><a href = {`https://www.metabolomicsworkbench.org/databases/refmet/refmet_details.php?REFMET_NAME=${val.REFMET_NAME}`} target = "_blank">{val.REFMET_NAME}</a></td>
              <td>
                {study_id_arr.map((study_id: string, i: number) =>
                  <>
                  {i > 0 ? <span>, </span> : null}
                  <a key={i} href = {`https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=${study_id}`} target = "_blank">{study_id}</a>
                  </>
                )}
              </td>
              </>
            </tr>
          )
        })}
        </table>
      </div>
    
    )
  })
  .build()
