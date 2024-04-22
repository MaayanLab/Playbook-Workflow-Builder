import { MetaNode } from '@/spec/metanode'
import { GeneSet } from '@/components/core/set'
import { MetgeneMetaboliteTable } from '../metgene_metabolite_table'
import { MetaboliteSet } from '@/components/core/set'
import { metgene_icon } from '@/icons'
import { MetGeneSummary } from '@/components/MW/metgene_summary'
import * as array from '@/utils/array'


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
  .output(MetgeneMetaboliteTable)  // Changed the codec to MGMetTableInfoC
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    const species_id = "hsa"
    const geneID_type = "SYMBOL_OR_ALIAS"
    const gene_ID = props.inputs.summary.gene

    // Original with vtf = "json"
    //const vtf = "json"; // vtf
    //const req = await fetch(`https://bdcw.org/MetGENE/rest/metabolites/species/${species_id}/GeneIDType/${geneID_type}/GeneInfoStr/${gene_ID}/anatomy/NA/disease/NA/phenotype/NA/viewType/${vtf}`)
    //const res = await req.json()
    //return res

    const vtf = "jsonfile"; // vtf
    const req = await fetch(`https://bdcw.org/MetGENE/rest/metabolites/species/${species_id}/GeneIDType/${geneID_type}/GeneInfoStr/${gene_ID}/anatomy/NA/disease/NA/phenotype/NA/viewType/${vtf}`)
    const res = await req.json();
    const jsonfileURL = res[0].FileURL;
    const datareq = await fetch(jsonfileURL);
    const data = await datareq.json();
    const retobj = {"jsonfile": {"FileDescription": "MetGENE metabolites json file", "FileURL": jsonfileURL},
          "contents": data}
    return retobj;
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
    const debug = 0;
    const species_id = "hsa"
    const geneID_type = "SYMBOL_OR_ALIAS"
    const gene_ID = props.inputs.geneset.set.join(",");

    const mw_dcc_cfde_server_url_dev = "https://sc-cfdewebdev.sdsc.edu";
    const mw_dcc_cfde_server_url_prod = "https://bdcw.org";
    const server_url = mw_dcc_cfde_server_url_prod; // mw_dcc_cfde_server_url_dev; 

    // Original with vtf = "json"
    //const vtf = "json"
    //const req = await fetch(`https://bdcw.org/MetGENE/rest/metabolites/species/${species_id}/GeneIDType/${geneID_type}/GeneInfoStr/${gene_ID}/anatomy/NA/disease/NA/phenotype/NA/viewType/${vtf}`)
    //const res = await req.json()
    //return  res

    const vtf = "jsonfile"; // vtf
    const requrl = `${server_url}/MetGENE/rest/metabolites/species/${species_id}/GeneIDType/${geneID_type}/GeneInfoStr/${gene_ID}/anatomy/NA/disease/NA/phenotype/NA/viewType/${vtf}`;
    if(debug > 0) {console.log(requrl);}
    const req = await fetch(requrl)
    if(debug > 0) {console.log(req);}
    const res = await req.json();
    if(debug > 0) {console.log(res);}
    const jsonfileURL = res[0].FileURL;
    if(debug > 0) {console.log(jsonfileURL);}
    const datareq = await fetch(jsonfileURL);
    const data = await datareq.json();
    const retobj = {"jsonfile": {"FileDescription": "MetGENE metabolites json file", "FileURL": jsonfileURL},
          "contents": data}
    return retobj;
  })
  .story(props =>
    `The gene set was then searched in the Metabolomics Workbench [Metabolomics Workbench, \\ref{https://www.metabolomicsworkbench.org/}] to identify associated metabolites.`
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
    // Mano: 2023/08/03: added contents since codec of MetgeneMetaboliteTable changed 
    const arr2 = props.inputs.mt.contents; 

    var MetArr = []; // Mano: Added 2023/06/28

    for (let i=0; i<arr2.length; i++) {
      for (let j=0; j<arr2[i].length; j++) {
        let metname_or_id = arr2[i][j].REFMET_NAME;
        if(! (metname_or_id  == "")){
          MetArr.push(metname_or_id ); // Added 2023/01/25
        }
      }
    }

    // keep only unique metabolites
    MetArr = array.unique(MetArr);
    return {"description": "", "set": MetArr} ; // Added 2023/01/25
  })
  .story(props =>
    `Then, MetaboliteSet (REFMET names) is extracted from the table of MetGENE Metabolites for various genes.`
  )
  .build()

// Mano: 2023/08/01: While keeping the original, writing a generalized/overloaded one
// Mano: Added 2023/06/28
export const MGMetTable2MetSet_T = [
  { id: 'REFMET_NAME' },
  { id: 'KEGG_COMPOUND_ID' },
].map(({id}) => MetaNode(`MGMetTable2MetSet_[${id}]`)
// Human readble descriptors about this node should go here
.meta({
  label: `Given the MetGENE Metabolites Table, generate a MetaboliteSet in terms of ${id}`,
  description: `Given the MetGENE Metabolites Table, generate a MetaboliteSet in terms of ${id}`,
})
// This should be a mapping from argument name to argument type
//  the types are previously defined Meta Node Data Types
.inputs({ mt: MetgeneMetaboliteTable })
// This should be a single Meta Node Data Type
.output(MetaboliteSet)
// The resolve function uses the inputs and returns output
//  both in the shape prescribed by the data type codecs
.resolve(async (props) => {
  const arr2 = props.inputs.mt.contents;

  const MetArr = []; // Mano: Added 2023/06/28

  for (let i=0; i<arr2.length; i++) {
    for (let j=0; j<arr2[i].length; j++) {
      //let metname_or_id = arr2[i][j].REFMET_NAME;
      
      let metname_or_id = "";
      if (id == 'REFMET_NAME') { 
        metname_or_id = arr2[i][j].REFMET_NAME;
      } else if (id == 'KEGG_COMPOUND_ID') {
        metname_or_id = arr2[i][j].KEGG_COMPOUND_ID;
      }
      if(! (metname_or_id  == "")){
        MetArr.push(metname_or_id ); // Added 2023/01/25
      }
    }
  }

  // keep only unique metabolites
  return {"description": `Metabolites (${id})`, "set": array.unique(MetArr)} ; // Added 2023/01/25
})
.story(props =>
  `Then, MetaboliteSet (${id}) is extracted from the table of MetGENE Metabolites for various genes.`
)
.build()
) // matches map
// The two resolver names will be: MGMetTable2MetSet_REFMET_NAME & MGMetTable2MetSet_KEGG_COMPOUND_ID
