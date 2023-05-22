import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { additional_info_icon, gene_icon } from '@/icons'
//import { VALID_LOADERS } from 'next/dist/shared/lib/image-config'

// How the schema validation works: https://codex.so/zod-validation-en
export const MyGeneIDOutputC = z.object({
  SYMBOL_OR_ALIAS: z.string().optional(),
  MATCH_SOURCE: z.string().optional(),
  ALIAS: z.string().optional(),
  ENSEMBL: z.string(),
  ENTREZID: z.string(),
  GENENAME: z.string(),
  REFSEQ: z.string().optional(),
  SYMBOL: z.string(),
  UNIPROT: z.string().optional(),
  KEGG: z.string().optional(),
  MARRVEL: z.string().optional()
})

export type MyGeneIDOutput = z.infer<typeof MyGeneIDOutputC>

// important ref: https://rsinohara.github.io/json-to-zod-react/
export const MyGeneIDOutputArrayC = z.array(
  MyGeneIDOutputC
)

export type MyGeneIDOutputArray = z.infer<typeof MyGeneIDOutputArrayC>

export function uniqJsonSubset(data:MyGeneIDOutputArray) {
  //let datastr = JSON.stringify(data);  datastr = datastr.replace("^[","");  datastr = datastr.replace("]$","");  const dataobj = JSON.parse(datastr);

  //const [dataobj] = data; // did not work
  //const dataobj = data; //const dataobj = [data][0];

  let uniqENTREZID : string[] = Array.from(new Set( data.map(a => a.ENTREZID)));
  let uniqSYMBOL : string[] = Array.from(new Set( data.map(a => a.SYMBOL)));
  let uniqGENENAME : string[] = Array.from(new Set( data.map(a => a.GENENAME)));

  let uniqdataobj = uniqENTREZID.map((id, idx) => {
   return {
    ENTREZID : uniqENTREZID[idx],
    SYMBOL : uniqSYMBOL[idx],
    GENENAME : uniqGENENAME[idx],
  };
  });
  return uniqdataobj;
}

// A unique name for your data type is used here
export const ConvertedGeneID = MetaNode('ConvertedGeneID')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'Gene IDs table',
    description: 'Gene IDs table',
    icon: [gene_icon, additional_info_icon],
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using zod, a compile-time and runtime type-safe codec can be constructed
  .codec(MyGeneIDOutputArrayC)
  // react component rendering your data goes here
  .view(data => {
    //return (...data) => {
    //window.alert("Type of data is:" + typeof(data));
    //const dataobj = [data][0];
    const dataobj = data; //const dataobj = [data][0];
    //window.alert(dataobj[0]);

    //for (let i=0; i<arr.length; i++) {document.writeln(arr[i]);}
    // extract all SYMBOL and ENTREZID: https://stackoverflow.com/questions/19590865/from-an-array-of-objects-extract-value-of-a-property-as-array
    // var, let, const: https://www.tutorialsteacher.com/typescript/typescript-variable
    // not used but keeping for knowing the syntax
    //let ENTREZID : string[] = dataobj.map(a => a.ENTREZID);
    //let SYMBOL : string[] = dataobj.map(a => a.SYMBOL);
    //let GENENAME : string[] = dataobj.map(a => a.GENENAME);

    // const picked = (({ a, c }) => ({ a, c }))(object);
    //let dataobj1 = Array.from(new Set(dataobj.map(a => ({"ENTREZID": a.ENTREZID,
    //                  "SYMBOL" : a.SYMBOL, "GENENAME": a.GENENAME}))) ];

    let uniqdataobj = uniqJsonSubset(data);

    return(
          <div><table>
    <tr>
    <th>SYMBOL</th>
    <th>NCBI ENTREZ URL</th>
    <th>GENE NAME</th>
    </tr>
{uniqdataobj.map((val, key) => {
      return (
      <tr key={key}>
        <td>{val.SYMBOL}</td>
        <td><a href = {`https://www.ncbi.nlm.nih.gov/gene/${val.ENTREZID}`}>{val.SYMBOL}</a></td>
        <td>{val.GENENAME}</td>
      </tr>
    )
    })}
    </table></div>
    )
  })
  .build()
