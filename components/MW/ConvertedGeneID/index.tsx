import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { VALID_LOADERS } from 'next/dist/shared/lib/image-config'

export function uniqJsonSubset(data) {
  let datastr = JSON.stringify(data);
  datastr = datastr.replace("^[","");
  datastr = datastr.replace("]$","");

  const dataobj = JSON.parse(datastr);
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
    let datastr = JSON.stringify(data);
    datastr = datastr.replace("^[","");
    datastr = datastr.replace("]$","");

    const dataobj = JSON.parse(datastr);
    //for (let i=0; i<arr.length; i++) {document.writeln(arr[i]);}
    //document.writeln(dataobj.length);
    //document.writeln(dataobj[0]["ENTREZID"]);
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

    let str = '';
    let stru = '';
    //var str1 : string[];
    for (var i = 0; i < ENTREZID.length; i++) {
        //str = str + "<a href=https://www.ncbi.nlm.nih.gov/gene/" + ENTREZID[i] + ">" + SYMBOL[i] + "</a> " + GENENAME[i] + "<br />";
        str = str + `<a href="https://www.ncbi.nlm.nih.gov/gene/${ENTREZID[i]}">${SYMBOL[i]}</a> ${GENENAME[i]}`;
    }
    for (var i = 0; i < uniqENTREZID.length; i++) {
      stru = stru + `<a href='https://www.ncbi.nlm.nih.gov/gene/${uniqENTREZID[i]}'>${uniqSYMBOL[i]}</a> ${uniqGENENAME[i]}`;
      //str1[i] = `<a href="https://www.ncbi.nlm.nih.gov/gene/${ENTREZID[i]}">${SYMBOL[i]}</a> ${GENENAME[i]}`;
    }
    return(
      //JSON.stringify(uniqdataobj)
          <div><table>
    <tr>
    <th>SYMBOL</th>
    <th>NCBI ENTREZ URL</th>
    <th>GENE NAME</th>
    </tr>
{uniqdataobj.map((val, key) => {
      return (//<div>junk</div>
      <tr key={key}>
        <td>{val.SYMBOL}</td>
        <td><a href = {`https://www.ncbi.nlm.nih.gov/gene/${val.ENTREZID}`}>{val.SYMBOL}</a></td>
        <td>{val.GENENAME}</td>
      </tr>
//<div>{JSON.stringify(data)}</div>
      //<div><h6>{str1.join(" ; ")}</h6></div>
      //{stru.map((val, key) => { 
      //    return (
      //      <div><h6><table><tr><td>{stru}</td></tr></table></h6></div> // ask Daniel why this str is not rendering as a web page/url
      //    )
    //})}
      //junk
      //<div>        <a href={`https://www.ncbi.nlm.nih.gov/gene/${geneinfo.entrezgene}`}>{geneinfo.symbol}</a> {geneinfo.name}      </div>
    ) //; };
    })}
    </table></div>
    )
  })
  .build()
