import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'

// A unique name for your data type is used here
export const MetGeneStudyTable = MetaNode.createData('MetGeneStudyTable')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGene Studies Table',
    description: 'Studies table corresponding to gene',
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using zod, a compile-time and runtime type-safe codec can be constructed
  .codec(z.any())
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
          {data.map((val, key) => {
          var all_study_ids = val.STUDY_ID
                
          var study_id_arr = all_study_ids.split(", ")
          
          function wrapURL(study_id) {
            return "https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID="+study_id  
          }
          var study_id_urls = study_id_arr.map(wrapURL)
          var urls = []
          for (var i = 0; i < study_id_arr.length; i++){
            var study_id = study_id_arr[i]+" "
            var study_url = study_id_urls[i]
            
            urls.push(<a href = {study_url} target = "_blank">{study_id}</a>)
            
          }
          return (
            <tr key={key}>
              <td>{val.KEGG_COMPOUND_ID}</td>
              <td><a href = {`https://www.metabolomicsworkbench.org/databases/refmet/refmet_details.php?REFMET_NAME=${val.REFMET_NAME}`} target = "_blank">{val.REFMET_NAME}</a></td>
              <td>{urls}</td>
            </tr>
          )
        })}
        </table>
      </div>
    )
  })
  .build()
