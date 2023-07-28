import { MetaNode } from '@/spec/metanode'
import { GeneTerm } from '@/components/core/input/term'
import { GeneSet } from '@/components/core/input/set'
import { MetgeneMetaboliteTable } from '../metgene_metabolite_table'
import { MetGeneMetObj } from '../metgene_metabolite_table'
import { MetGeneMetObjArray } from '../metgene_metabolite_table'
import { MetaboliteSet } from '@/components/core/input/set'
import { metgene_icon } from '@/icons'
import { MetGeneSummary } from '../metgene_summary'


// A unique name for your resolver is used here
export const MetgeneMetabolites = MetaNode('MetgeneMetabolites')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGENE Metabolites',
    description: 'Extract Metabolomics metabolites for the gene from MetGENE',
    icon: [metgene_icon],
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ summary: MetGeneSummary })
  // This should be a single Meta Node Data Type
  .output(MetgeneMetaboliteTable)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    const species_id = "hsa"
    const geneID_type = "SYMBOL_OR_ALIAS"
    const gene_ID = props.inputs.summary.gene
    const vtf = "json"
    const req = await fetch(`https://bdcw.org/MetGENE/rest/metabolites/species/${species_id}/GeneIDType/${geneID_type}/GeneInfoStr/${gene_ID}/anatomy/NA/disease/NA/phenotype/NA/viewType/${vtf}`)
    const res = await req.json()


    //return props.inputs.input
    return  res
  })
  .story(props =>
    `${props.inputs ? props.inputs.summary.gene : 'The gene'} was then searched in the Metabolomics Workbench [\\ref{The Metabolomics Workbench, https://www.metabolomicsworkbench.org/}] to identify associated metabolites.`
  )
  .build()

export const MetgeneMetabolitesGeneSet = MetaNode('MetgeneMetabolitesGeneSet')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'MetGENE Metabolites with Gene Set',
    description: 'Compute the MetGENE metabolites for a Gene Set',
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
    const geneID_type = "SYMBOL_OR_ALIAS"
    const gene_ID = props.inputs.geneset.set.join(",");
    const vtf = "json"
    const req = await fetch(`https://bdcw.org/MetGENE/rest/metabolites/species/${species_id}/GeneIDType/${geneID_type}/GeneInfoStr/${gene_ID}/anatomy/NA/disease/NA/phenotype/NA/viewType/${vtf}`)
    const res = await req.json()


    //return props.inputs.input
    return  res
  })
  .story(props =>
    `The gene set was then searched in the Metabolomics Workbench [REF] to identify associated metabolites.`
  )
  .build()

// Mano: Added 2023/06/28
  export const MGMetTable2MetSet = MetaNode('MGMetTable2MetSet')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'Given the MetGENE Metabolites Table, generate a MetaboliteSet',
    description: 'Given the MetGENE Metabolites Table, generate a MetaboliteSet',
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ mt: MetgeneMetaboliteTable })
  // This should be a single Meta Node Data Type
  .output(MetaboliteSet)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    const arr2 = props.inputs.mt;

    var MetArr = []; // Mano: Added 2023/06/28

    for (let i=0; i<arr2.length; i++) {
      for (let j=0; j<arr2[i].length; j++) {
        if(! (arr2[i][j].REFMET_NAME == "")){
          MetArr.push(arr2[i][j].REFMET_NAME); // Added 2023/01/25
        }
      }
    }

    // keep only unique metabolites
    MetArr = [...new Set( MetArr )];
    return {"description": "", "set": MetArr} ; // Added 2023/01/25
  })
  .story(props =>
    `Then, MetaboliteSet (REFMET names) is extracted from the table of MetGENE Metabolites for various genes.`
  )
  .build()
