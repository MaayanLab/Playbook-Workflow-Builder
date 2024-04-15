import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { additional_info_icon } from '@/icons';

// Mano: 2023/08/01: adding more fields and making some optional
const MetSummaryObjC = z.object({
  name: z.string().optional(),
  refmet_name: z.string().optional(),
  pubchem_cid: z.string().optional(),
  inchi_key: z.string().optional(),
  exactmass: z.string().optional(),
  formula: z.string().optional(),
  super_class: z.string().optional(),
  main_class: z.string().optional(),
  sub_class: z.string().optional(),
});

export type  MetSummaryObj = z.infer<typeof  MetSummaryObjC>

export const MetSummaryObjArrayC = z.array(
  MetSummaryObjC
)

export type MetSummaryObjArray = z.infer<typeof MetSummaryObjArrayC>


// A unique name for your data type is used here
export const MetaboliteSummary = MetaNode('MetaboliteSummary')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'Metabolite Summary',
    description: 'Metabolite Summary, rendered',
    icon: [additional_info_icon],
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using zod, a compile-time and runtime type-safe codec can be constructed
  .codec(MetSummaryObjArrayC)
  // react component rendering your data goes here
  .view(data => {
    //var data1 = [data];
    return (
      //<div>{JSON.stringify(data)}</div>
      <div className="prose max-w-none">
        <h2>Metabolite Summary</h2>
        <table>
          <th> Name </th>
          <th> RefMet Name </th>
          <th> PubChem CID</th>
          <th> Exact Mass</th>
          <th> Formula </th>
          <th> Super class </th>
          <th> Main class </th>
          <th> Sub class </th>
          <th> InChI</th>
        {data.map((val: MetSummaryObj, i: number) => {
            //window.alert(val); // console.log(val) // //window.alert(dataobj[0]);
            return (
              <tr key={i}>
                <td>{val.name}</td>
                <td>{val.refmet_name}</td>
                <td>{val.pubchem_cid}</td>
                <td>{val.exactmass}</td>
                <td>{val.formula}</td>
                <td>{val.super_class}</td>
                <td>{val.sub_class}</td>
                <td>{val.main_class}</td>
                <td>{val.inchi_key}</td>
              </tr>
            )
        })}



        </table>
      </div>
    )
  })
  .build()
