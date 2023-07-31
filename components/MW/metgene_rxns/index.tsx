import { MetaNode } from '@/spec/metanode'
import { GeneSet } from '@/components/core/input/set'
import { MetGeneRxnTable } from '@/components/MW/metgene_rxn_table'
import { metgene_icon } from '@/icons'
import { MetGeneSummary } from '@/components/MW/metgene_summary'

// A unique name for your resolver is used here
export const MetGeneRxns = MetaNode('MetGeneRxns')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGENE Reactions',
    description: 'Extract Metabolomics reactions for the gene from MetGENE',
    icon: [metgene_icon],
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ summary: MetGeneSummary })
  // This should be a single Meta Node Data Type
  .output(MetGeneRxnTable)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    const species_id = "hsa"
    const geneID_type = "SYMBOL_OR_ALIAS"

    const gene_ID = props.inputs.summary.gene
    const vtf = "json"
    const req = await fetch(`https://bdcw.org/MetGENE/rest/reactions/species/${species_id}/GeneIDType/${geneID_type}/GeneInfoStr/${gene_ID}/anatomy/NA/disease/NA/phenotype/NA/viewType/${vtf}`)
    const res = await req.json()


    //return props.inputs.input
    return  res

  })
  .story(props =>
    `${props.inputs ? props.inputs.summary.gene : 'The gene'} was then searched in the Metabolomics Workbench [\\ref{The Metabolomics Workbench, https://www.metabolomicsworkbench.org/}] to identify relevant reactions.`
  )
  .build()

export const MetGeneRxnsGeneSet = MetaNode('MetGeneRxnsSet')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGENE Reactions with Gene Set',
    description: 'Compute the MetGENE Reactions for a Gene Set',
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ geneset: GeneSet })
  // This should be a single Meta Node Data Type
  .output(MetGeneRxnTable)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    const species_id = "hsa"
    const geneID_type = "SYMBOL_OR_ALIAS"

    const gene_ID = props.inputs.geneset.set.join(",");
    const vtf = "json"
    const req = await fetch(`https://bdcw.org/MetGENE/rest/reactions/species/${species_id}/GeneIDType/${geneID_type}/GeneInfoStr/${gene_ID}/anatomy/NA/disease/NA/phenotype/NA/viewType/${vtf}`)
    const res = await req.json()


    //return props.inputs.input
    return  res

  })
  .story(props =>
    `The gene set was then searched in the Metabolomics Workbench [\\ref{The Metabolomics Workbench, https://www.metabolomicsworkbench.org/}] to identify relevant reactions.`
  )
  .build()
