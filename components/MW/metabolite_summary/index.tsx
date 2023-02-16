import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'

const MetSummaryObjC = z.object({
  name: z.string(),
  pubchem_cid: z.string(),
  inchi_key: z.string(),
  exactmass: z.string(),
  formula: z.string(),
  super_class: z.string(),
  main_class: z.string(),
  sub_class: z.string(),
});

export type  MetSummaryObj = z.infer<typeof  MetSummaryObjC>

// A unique name for your data type is used here
export const MetaboliteSummary = MetaNode('MetaboliteSummary')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'Metabolite Summary',
    description: 'Metabolite Summary, rendered',
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using zod, a compile-time and runtime type-safe codec can be constructed
  .codec(MetSummaryObjC)
  // react component rendering your data goes here
  .view(data => {
    var data1 = [data];
    return (
      //<div>{JSON.stringify(data)}</div>
      <div>
        <h2>Metabolite Summary</h2>
        <table>
          <th> Name </th>
          <th> PubChem ID</th>
          <th> Exact Mass</th>
          <th> Fromula </th>
          <th> Super class </th>
          <th> Main class </th>
          <th> Sub class </th>
          <th> InChI</th>
        {data1.map((val: MetSummaryObj, key: number) => {
          return (
            <tr key={key}>
              <td>{val.name}</td>
              <td>{val.pubchem_cid}</td>
              <td>{val.exactmass}</td>
              <td>{val.formula}</td>
              <td>{val.super_class}</td>
              <td>{val.sub_class}</td>
              <td>{val.main_class}</td>
              <td>{val.inchi_key}</td>
            </tr>
          )})
        }
        </table>
      </div>
    )
  })
  .build()
