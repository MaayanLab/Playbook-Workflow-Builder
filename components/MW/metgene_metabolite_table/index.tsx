import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { FileDesURLC } from '../MetENP_on_MetSet'
import { additional_info_icon, drug_icon } from '@/icons';

const MetGeneMetObjC = z.object({
  Gene: z.string(),
  KEGG_COMPOUND_ID: z.string(),
  REFMET_NAME: z.string(),
  KEGG_REACTION_ID: z.string(),
  METSTAT_LINK: z.string(),
});

export type  MetGeneMetObj = z.infer<typeof  MetGeneMetObjC>

// important ref: https://rsinohara.github.io/json-to-zod-react/
// array of objects: [{}]
export const MetGeneMetObjArrayC = z.array(
  MetGeneMetObjC
)

export type MetGeneMetObjArray = z.infer<typeof MetGeneMetObjArrayC>

// array of array of objects: [[{}]]
export const MetGeneMetObjArray2C = z.array(
  MetGeneMetObjArrayC
)

export type MetGeneMetObjArray2 = z.infer<typeof MetGeneMetObjArray2C>

// Mano: 2023/08/01: To allow getting just the file URL (the file will have the actual 2D-array json content)
// Schema for output from application of MetENP on a list of metabolites (MetaboliteSet)
export const MGMetTableInfoC = z.object({
  jsonfile: FileDesURLC,
  contents: MetGeneMetObjArray2C
})
export type MGMetTableInfo = z.infer<typeof MGMetTableInfoC>

// A unique name for your data type is used here
export const MetgeneMetaboliteTable = MetaNode('MetgeneMetaboliteTable')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGENE metabolite table',
    description: 'MetGENE metabolite table',
    icon: [drug_icon, additional_info_icon],
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using zod, a compile-time and runtime type-safe codec can be constructed
  //.codec(MetGeneMetObjArray2C) // Mano: 2023/08/03: Original
  .codec(MGMetTableInfoC)
  // react component rendering your data goes here
  .view(data => {
    const heading1 = "Gene"
    const heading2 = "KEGG_COMPOUND_ID"
    const heading3 = "REFMET_NAME"
    const heading4 = "KEGG_REACTION_ID"
    const heading5 = "METSTAT_LINK"

    const MaxRows2Show = 50;
    const ng = data.contents.length;
    const arrlen:number[] = data.contents.map(x=>x.length);
    let arrlen_cumsum: number[] = [];  
    arrlen.reduce( (prev, curr,i) =>  arrlen_cumsum[i] = prev + curr , 0 );
    //console.log(arrlen); console.log(arrlen_cumsum);
    //window.alert(arrlen); window.alert(arrlen_cumsum);

    let isLessThanMax = (x: number) => x < MaxRows2Show;
    //didn't work: let ng2show = arrlen_cumsum.findLastIndex(isLessThanMax); if(ng2show < 0) ng2show = 0; // Math.min(MaxRows2Show, ng);
    let ng2show = 0; for(let ii = 0; ii < arrlen_cumsum.length; ii++){ ng2show = ii; if(arrlen_cumsum[ii] > MaxRows2Show) {break;}}
    // Show up to the point, total numrows so far is about less than MaxRows2Show, OK to show it for next gene but no more

    ng2show = ng2show + 1; // actual number rather than index
    let ExtraText_if_more_rows = "";
    if(ng2show < ng) ExtraText_if_more_rows = ` The tables for the first ${ng2show} genes are shown below.`;

    // Below, data.contents is array of array of objects; trying to show at most 50 rows
    return (
      <div>
        <h2>MetGENE metabolites</h2>
        <span style={{ color: "#0000F0" }}>
          <h3>Full table as a json file is available at <a href={`${data.jsonfile.FileURL}`} target="__blank">{data.jsonfile.FileDescription}</a>. {ExtraText_if_more_rows}</h3>
        </span>
        {data.contents.slice(0, ng2show).map((arrayVal:MetGeneMetObjArray, index:number) => (
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

                {arrayVal.map((val:MetGeneMetObj, j:number) => {
                  const all_rxn_ids = val.KEGG_REACTION_ID
                  const cpd_id = val.KEGG_COMPOUND_ID
                  const rxn_id_arr = all_rxn_ids.split(", ")


                    return (
                      <tr key={j}>
                        <td>{val.Gene}</td>
                        <td><a href = {`https://www.kegg.jp/entry/${val.KEGG_COMPOUND_ID}`} target = "__blank">{val.KEGG_COMPOUND_ID}</a></td>
                        <td><a href = {`https://www.metabolomicsworkbench.org/databases/refmet/refmet_details.php?REFMET_NAME=${val.REFMET_NAME}`} target = "_blank">{val.REFMET_NAME}</a></td>
                        <td>
                          {rxn_id_arr.map((rxn_id:string, i:number) =>
                            <>
                            {i > 0 ? <span>, </span> : null}
                            <a key={i} href = {`https://www.genome.jp/entry/rn:${rxn_id}`} target = "_blank">{rxn_id}</a>
                            </>
                          )}
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
