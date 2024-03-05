import React from "react";
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { glygen_icon } from '@/icons'
import { GlycanTerm } from '@/components/core/input/term'
import { GlycanSet } from '@/components/core/input/set'
import { GlycanResponse, GlyGenGlycanSetResponse } from "./data_models";
import { GlycanClassification, GlycanCrossRef } from "./sup_components";
import { glycan_set_search_query } from "./sup_functions";

// --- Data Metanodes --- //

// Single glycan data metanode 
export const GlycanViewResponseNode = MetaNode('GlycanViewResponse')
  .meta({
    label: 'Glycan information',
    description: 'Glycan information from GlyGen'
  })
  .codec(GlycanResponse)
  .view(data => {
    const glyGenLink = `http://www.glygen.org/glycan/${data.glytoucan.glytoucan_ac}`;

    return (
      <div className="prose">
        <div>GlyTouCan Accession: 
          <b>
            <a href={glyGenLink} target='_blank' rel='noopener nonreferrer' style={{color: 'blue'}}> <u style={{color: 'blue'}}>{data.glytoucan.glytoucan_ac}</u></a>
          </b>
        </div>
        <div>Monoisotopic Mass: <b>{data.mass} Da</b></div>
        <div>Monoisotopic Mass-pMe (Da): <b>{data.mass_pme} Da</b></div>
        <div>
          Glycan Type / Glycan Subtype: <b><GlycanClassification classification={data.classification}/></b>
        </div>
        <div>
          <GlycanCrossRef crossref={data.crossref}/>
        </div>
        <div>
          Glycan Image: 
          <img src={`https://api.glygen.org/glycan/image/${data.glytoucan.glytoucan_ac}/`} alt='Glycan Image'/>
        </div>
      </div>
    )
  })
  .build()

// Glycan set data metanode
export const GlyGenGlycanSetResponseNode = MetaNode('GlyGenGlycanSetResponse')
  .meta({
    label: 'GlyGen Glycans',
    description: 'Protein product records in GlyGen',
    icon: [glygen_icon],
  })
  .codec(GlyGenGlycanSetResponse)
  .view( data => {
    return (
      <div>
          <table>
            <thead>
              <tr>
                <th>Glycan ID</th>
                <th>Glycan Image</th>
                <th>Hit Score</th>
                <th>Monoisotopic Mass</th>
                <th>Monoisotopic Mass-pMe (Da)</th>
                <th>No of Sugars</th>
                <th>No of Glycoproteins</th>
              </tr>
            </thead>
            <tbody>
              {data.map((entry, index) => (
                <tr key = {index}>
                  <td>{entry.glytoucan.glytoucan_ac}</td>
                  <td><img src={`https://api.glygen.org/glycan/image/${entry.glytoucan.glytoucan_ac}/`} alt='Glycan Image'/></td>
                  <td>{entry.hit_score}</td>
                  <td>{`${entry.mass}`}</td>
                  <td>{entry.mass_pme}</td>
                  <td>{entry.sugar_count}</td>
                  <td>{entry.glycoprotein_count}</td>
                  <td>{entry.associated_enzymes}</td>
                </tr>
              ))}
            </tbody>
          </table>
      </div>
    )
  })
  .build()


// --- Process Metanodes --- //

// Single glycan process metanode
export const GlycanInformation = MetaNode('GlycanInformation')
  .meta({
    label: 'Search GlyGen by GlyTouCan Accession',
    description: 'Search for Glycan information',
    // icon: []
    pagerank: 2
  })
  .inputs({ glycan: GlycanTerm })
  .output(GlycanViewResponseNode)
  .resolve(async (props) => {
    // get glycan data 
    const detail_response = await fetch(`https://api.glygen.org/glycan/detail/${props.inputs.glycan}`, {
      method: 'POST',
      headers: {
        accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ glytoucan_ac: props.inputs.glycan })
    })
    const glycan_data = await detail_response.json();
    return glycan_data
  })
  .story(props => 
    `The GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a information about ${props.inputs ? props.inputs.glycan : 'the glycan'}.`
  )
  .build()

// Glycan set process metanode
export const GlyGenGlycanSet = MetaNode('GGGS')
  .meta({
    label: 'Search GlyGen by GlyTouCan Accession for Glycans',
    description: 'Find glycan records in GlyGen.',
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ glycan_glytoucan_acc_set: GlycanSet })
  .output(GlyGenGlycanSetResponseNode)
  .resolve(async (props) => {
    console.log("===> Got glycan set input: ", props.inputs.glycan_glytoucan_acc_set.set);
    const glycan_response = await glycan_set_search_query(props.inputs.glycan_glytoucan_acc_set.set);
    return glycan_response;
  })
  .story(props =>
    // TODO: re-write story sentence to make sense with protein term input (previous gene value removed to prevent `npm run build` error)
    `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of glycans that originate from.`
  )
  .build()

