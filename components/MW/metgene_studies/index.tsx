import { MetaNode } from '@/spec/metanode'
import { GeneSet } from '@/components/core/set'
import { MetGeneStudyTable } from '@/components/MW/metgene_study_table'
import { metgene_icon, additional_info_icon, gene_icon, set_icon } from '@/icons'
import { MetGeneSummary } from '@/components/MW/metgene_summary'
import Citable from '@/utils/citations'

// A unique name for your resolver is used here
export const MetGeneStudies = MetaNode('MetGeneStudies')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGENE Studies',
    description: 'Extract Metabolomics studies for the gene from MetGENE',
    icon: [metgene_icon],
  })


  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ summary: MetGeneSummary })
  // This should be a single Meta Node Data Type
  .output(MetGeneStudyTable)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    const species_id = "hsa"
    const geneID_type = "SYMBOL_OR_ALIAS"

    const gene_ID = props.inputs.summary.gene
    const vtf = "json"
    const req = await fetch(`https://bdcw.org/MetGENE/rest/studies/species/${species_id}/GeneIDType/${geneID_type}/GeneInfoStr/${gene_ID}/anatomy/NA/disease/NA/phenotype/NA/viewType/${vtf}`)
    const res = await req.json()



    return  res
  })
  .story(props => ({
    abstract: Citable.text`${props.inputs?.summary?.gene ? props.inputs.summary.gene : 'The gene'} was then searched in the Metabolomics Workbench [${Citable.cite('The Metabolomics Workbench, https://www.metabolomicsworkbench.org/')}] to identify relevant studies related to the gene.`
  }))
  .build()

export const MetGeneStudiesGeneSet = MetaNode('MetGeneStudiesGeneSet')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGENE Studies with gene set',
    description: 'Compute the MetGENE studies function for a gene set',
    icon: [gene_icon, set_icon, additional_info_icon],
  })


  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ geneset: GeneSet })
  // This should be a single Meta Node Data Type
  .output(MetGeneStudyTable)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    const species_id = "hsa"
    const geneID_type = "SYMBOL_OR_ALIAS"

    const gene_ID = props.inputs.geneset.set.join(",");
    const vtf = "json"
    const req = await fetch(`https://bdcw.org/MetGENE/rest/studies/species/${species_id}/GeneIDType/${geneID_type}/GeneInfoStr/${gene_ID}/anatomy/NA/disease/NA/phenotype/NA/viewType/${vtf}`)
    const res = await req.json()



    return  res
  })
  .story(props => ({
    abstract: Citable.text`The gene set was then searched in the Metabolomics Workbench [${Citable.cite('The Metabolomics Workbench, https://www.metabolomicsworkbench.org/')}] to identify relevant studies related to the genes.`
  }))
  .build()