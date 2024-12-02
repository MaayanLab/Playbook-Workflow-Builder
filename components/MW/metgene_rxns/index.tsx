import { MetaNode } from '@/spec/metanode'
import { GeneSet } from '@/components/core/set'
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
    external: true,
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
  .story(props => ({
    abstract: `${props.inputs?.summary?.gene ? props.inputs.summary.gene : 'The gene'} was then searched in the Metabolomics Workbench\\ref{The Metabolomics Workbench, https://www.metabolomicsworkbench.org/} to identify relevant reactions.`,
    introduction: `MetGENE is a information retrieval tool that connects a gene or a set of genes to metabolomic studies in the Metabolomic Workbench. It uses a knowledge based approach where the gene is connected to pathways it regulates, followed by reactions within the pathways and metabolites particicpating in the reactions. The metabolites are connected to studies in Metabolomics Workbench.`,
    methods: `Given a gene, MetGENE provides REST API to extract information regarding all the reactions in the pathways where the gene participates e.g.  for human species (hsa), with anatomy blood and disease diabetes,\\ref{https://bdcw.org/MetGENE/rest/reactions/species/hsa/GeneIDType/SYMBOL/GeneInfoStr/HK1/anatomy/blood/disease/diabetes/phenotype/NA/viewType/json}, returns the KEGG Reaction IDS, KEGG Reaction names as well as KEGG Reaction equations for the corresponding gene in the form of a JSON output. `,
    legend: `The list of results are displayed as a table comprising of KEGG Reaction IDs linked to corresponding reaction descriptions, KEGG Reaction names and KEGG Reaction Equations from the KEGG database. The reaction IDs are hyperlinked to the respective reaction information page on the KEGG web site.`,
  }))
  .build()

export const MetGeneRxnsGeneSet = MetaNode('MetGeneRxnsSet')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGENE Reactions with Gene Set',
    description: 'Compute the MetGENE Reactions for a Gene Set',
    external: true,
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
  .story(props => ({
    abstract: `The gene set was then searched in the Metabolomics Workbench\\ref{The Metabolomics Workbench, https://www.metabolomicsworkbench.org/} to identify relevant reactions.`,
    introduction: `MetGENE is a information retrieval tool that connects a gene or a set of genes to metabolomic studies in the Metabolomic Workbench. It uses a knowledge based approach where the gene is connected to pathways it regulates, followed by reactions within the pathways and metabolites particicpating in the reactions. The metabolites are connected to studies in Metabolomics Workbench.`,
    methods: `Given a list of genes, MetGENE provides REST API to extract information regarding all the reactions in the pathways where the genes participatese.g.  for a given species (hsa), with anatomy blood and disease diabetes,\\ref{https://bdcw.org/MetGENE/rest/metabolites/species/hsa/GeneIDType/SYMBOL/GeneInfoStr/HK1,RPE/anatomy/NA/disease/Diabetes/phenotype/NA/viewType/json}, returns the KEGG Reaction IDS, KEGG Reaction names as well as KEGG Reaction equations for the corresponding gene in the form of a JSON output.`,
    legend: `For each gene in the geneset, the list of results are displayed in the table format, comprising of KEGG reaction IDs corresponding to the reactions, KEGG Reaction name as well as the KEGG Reaction Equation from the KEGG database. The reaction IDs are hyperlinked to the respective reaction information page on the KEGG web site.`,
  }))
  .build()
