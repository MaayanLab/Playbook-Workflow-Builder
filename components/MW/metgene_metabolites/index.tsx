import { MetaNode } from '@/spec/metanode'
import { GeneTerm } from '@/components/core/input/term'
import { GeneSet } from '@/components/core/input/set'
import { MetgeneMetaboliteTable } from '../metgene_metabolite_table'


// A unique name for your resolver is used here
export const MetgeneMetabolites = MetaNode.createProcess('MetgeneMetabolites')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGENE Metabolites',
    description: 'Compute the MetGENE metabolites',
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ gene: GeneTerm })
  // This should be a single Meta Node Data Type
  .output(MetgeneMetaboliteTable)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    const species_id = "hsa"
    const geneID_type = "SYMBOL"
    const gene_ID = props.inputs.gene
    const vtf = "json"
    const req = await fetch(`https://bdcw.org/MetGENE/rest/metabolites/species/${species_id}/GeneIDType/${geneID_type}/GeneInfoStr/${gene_ID}/anatomy/NA/disease/NA/phenotype/NA/viewType/${vtf}`)
    const res = await req.json()
    

    //return props.inputs.input
    return  res
  })
  .build()

  export const MetgeneMetabolitesGeneSet = MetaNode.createProcess('MetgeneMetabolitesGeneSet')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGENE Metabolites with GeneSet',
    description: 'Compute the MetGENE metabolites for a GeneSet',
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ geneset: GeneSet })
  // This should be a single Meta Node Data Type
  .output(MetgeneMetaboliteTable)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    const species_id = "hsa"
    const geneID_type = "SYMBOL"
    const gene_ID = props.inputs.geneset.join(",");
    const vtf = "json"
    const req = await fetch(`https://bdcw.org/MetGENE/rest/metabolites/species/${species_id}/GeneIDType/${geneID_type}/GeneInfoStr/${gene_ID}/anatomy/NA/disease/NA/phenotype/NA/viewType/${vtf}`)
    const res = await req.json()
    

    //return props.inputs.input
    return  res
  })
  .build()