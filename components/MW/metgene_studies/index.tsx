import { MetaNode } from '@/spec/metanode'
import { GeneTerm } from '@/components/core/input/term'

import { MetGeneStudyTable } from '../metgene_study_table'

// A unique name for your resolver is used here
export const MetGeneStudies = MetaNode('MetGeneStudies')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGENE Studies',
    description: 'Compute the MetGENE studies function',
  })


  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ gene: GeneTerm })
  // This should be a single Meta Node Data Type
  .output(MetGeneStudyTable)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    const species_id = "hsa"
    const geneID_type = "SYMBOL"
    
    const gene_ID = props.inputs.gene
    const vtf = "json"
    const req = await fetch(`https://bdcw.org/MetGENE/rest/studies/species/${species_id}/GeneIDType/${geneID_type}/GeneInfoStr/${gene_ID}/anatomy/NA/disease/NA/phenotype/NA/viewType/${vtf}`)
    const res = await req.json()
    

    
    return  res
  })
  .build()
