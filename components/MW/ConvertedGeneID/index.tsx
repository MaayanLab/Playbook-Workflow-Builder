import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { VALID_LOADERS } from 'next/dist/shared/lib/image-config'

export function uniqJsonSubset(data) {
  //let datastr = JSON.stringify(data);  datastr = datastr.replace("^[","");  datastr = datastr.replace("]$","");  const dataobj = JSON.parse(datastr);
  
  const [dataobj] = data;
  let uniqENTREZID : string[] = Array.from(new Set( dataobj.map(a => a.ENTREZID)));
  let uniqSYMBOL : string[] = Array.from(new Set( dataobj.map(a => a.SYMBOL)));
  let uniqGENENAME : string[] = Array.from(new Set( dataobj.map(a => a.GENENAME)));

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
export const ConvertedGeneID = MetaNode.createData('ConvertedGeneID')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'Json for converted gene IDs table',
    description: 'Json for converted gene IDs table',
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using zod, a compile-time and runtime type-safe codec can be constructed
  .codec(z.any())
  // react component rendering your data goes here
  .view(data => {
    //return (...data) => {
    const [dataobj] = data;

    //for (let i=0; i<arr.length; i++) {document.writeln(arr[i]);}
    // extract all SYMBOL and ENTREZID: https://stackoverflow.com/questions/19590865/from-an-array-of-objects-extract-value-of-a-property-as-array
    // var, let, const: https://www.tutorialsteacher.com/typescript/typescript-variable
    let ENTREZID : string[] = dataobj.map(a => a.ENTREZID);
    let SYMBOL : string[] = dataobj.map(a => a.SYMBOL);
    let GENENAME : string[] = dataobj.map(a => a.GENENAME);

    let uniqENTREZID : string[] = Array.from(new Set( dataobj.map(a => a.ENTREZID)));
    let uniqSYMBOL : string[] = Array.from(new Set( dataobj.map(a => a.SYMBOL)));
    let uniqGENENAME : string[] = Array.from(new Set( dataobj.map(a => a.GENENAME)));

    // const picked = (({ a, c }) => ({ a, c }))(object);
    //let dataobj1 = Array.from(new Set(dataobj.map(a => ({"ENTREZID": a.ENTREZID, 
    //                  "SYMBOL" : a.SYMBOL, "GENENAME": a.GENENAME}))) ];

    let uniqdataobj = uniqENTREZID.map((id, idx) => {
     return {
      ENTREZID : uniqENTREZID[idx], 
      SYMBOL : uniqSYMBOL[idx],
      GENENAME : uniqGENENAME[idx],
    };
    });

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
