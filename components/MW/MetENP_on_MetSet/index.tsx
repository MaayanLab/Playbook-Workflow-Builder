import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { MetaboliteSet } from '@/components/core/set'
import { metabolomicsworkbench_icon } from '@/icons'
import { z } from 'zod'
import { additional_info_icon} from '@/icons'

//import { VALID_LOADERS } from 'next/dist/shared/lib/image-config'


// How the schema validation works: https://codex.so/zod-validation-en

// Schema for output from application of MetENP on a list of metabolites (MetaboliteSet)
export const FileDesURLC = z.object({  FileDescription: z.string(), FileURL: z.string() })

export type FileDesURL = z.infer<typeof FileDesURLC>

export const MetENP_MetSet_outputC = z.object({
  RESTURL: z.string(), 
  res: z.array(FileDesURLC)
})

export type MetENP_MetSet_output = z.infer<typeof MetENP_MetSet_outputC>
  
// A unique name for your data type is used here
export const MetENP_MetSet = MetaNode('MetENP_MetSet')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'Output from MetENP on MetaboliteSet',
    description: 'Output from MetENP on MetaboliteSet',
    icon: [additional_info_icon],
    external: true,
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using zod, a compile-time and runtime type-safe codec can be constructed
  .codec(MetENP_MetSet_outputC)
  // react component rendering your data goes here
  .view(data => {
    const dataobj = data.res; //const dataobj = [data][0];

    return(
          <div className="prose max-w-none">
            <table>
    <tr>
    <th>No.</th>
    <th>Description/URL</th>
    </tr>
{dataobj.map((val, key) => {
      return (
      <tr key={key}>
        <td>{key}</td>
        <td><a href = {val.FileURL}>{val.FileDescription}</a></td>
      </tr>
    )
    })}
    </table></div>
    )
  })
  .build()

// A unique name for your resolver is used here
export const Call_MetENP_on_MetSet = MetaNode('Call_MetENP_on_MetSet')
// Human readble descriptors about this node should go here
.meta({
  label: 'Perform MetENP analysis on a list of metabolites (no numeric data)',
  description: 'Given a list of metabolites as MetaboliteSet, Perform MetENP analysis on the list',
  icon: [metabolomicsworkbench_icon],
  external: true,
})
// This should be a mapping from argument name to argument type
//  the types are previously defined Meta Node Data Types
.inputs({ mset: MetaboliteSet })
// This should be a single Meta Node Data Type
.output(MetENP_MetSet)
// The resolve function uses the inputs and returns output
//  both in the shape prescribed by the data type codecs
.resolve(async (props) => {
  // Example URL, triple encodes the list of metabolites to include in the URL
  // https://sc-cfdewebdev.sdsc.edu//MetENP/rest/metclass/sub_class/updown_fillcolor/red__green__blue/enrich_stats/HG/no/1/sps/hsa/padj/fdr/kegg_comp_path/FALSE/geneoption/TRUE/PPIopt/1_1000_400/location/0/metlistfname/Acetic%252520acid%25255CnGlucose%25255Cn11b%252520PGE2%25255Cn11%252520HDoHE%25255Cn11-HEPE%25255Cn12%25252C13%252520diHOME%25255CnFA%25252822%25253A1%252529%25255CnNADP%25252B%25255Cn3%252527-UMP
  // shorter list: https://sc-cfdewebdev.sdsc.edu//MetENP/rest/metclass/sub_class/updown_fillcolor/red__green__blue/enrich_stats/HG/no/1/sps/hsa/padj/fdr/kegg_comp_path/FALSE/geneoption/TRUE/PPIopt/1_1000_400/location/0/metlistfname/Glucose%25255Cn11b%252520PGE2%25255Cn11%252520HDoHE%25255Cn11-HEPE%25255Cn3%252527-UMP
  const species_txid = "9606"
  const species_id = "hsa"
  const string_api_url = "https://version-11-5.string-db.org/api"
  const output_format = "json"
  const method = "interaction_partners"

  //const metenp_api_url = "https://sc-cfdewebdev.sdsc.edu//MetENP/rest/";
  const metenp_api_url = "https://bdcw.org//MetENP/rest/";
  const metenp_basic_param1 = "metclass/sub_class/updown_fillcolor/red__green__blue/enrich_stats/HG/no/1/";
  const sps_str = "sps/" + species_id + "/";
  // can set geneoption to FALSE as TRUE option can take a long time (> 10 mins fr RPE-> PPI -> Metgene metabolites)
  const metenp_basic_param2 = "padj/fdr/kegg_comp_path/FALSE/geneoption/TRUE/"; // geneoption/TRUE
  // If don't want PPI, set the first parameter to 0 below
  const PPI_optstr = "PPIopt/1_1000_700/"; // "PPIopt/1_1000_700/";
  const location_optstr = "location/0/";
  const metlist_as_nlstr = props.inputs.mset.set.join("\\n");
  
  // to test, can use a list of metabolites:
  //const metlist_as_nlstr = ["Acetic acid", "Glucose", "Fructose"].join("\\n"); // use \\n not just \n
 
  // List of metabolites should be separated by \n (actually \\n in terms of code) and 
  // triple encoded since the server seems to be decoding it twice before it goes into the REST php code.
  const metlist_as_nlstr_enc3 = encodeURIComponent(encodeURIComponent(encodeURIComponent(metlist_as_nlstr)))
  const metlist_optstr = "metlistfname/" + metlist_as_nlstr_enc3;
  
  const request_url = [metenp_api_url, metenp_basic_param1, sps_str, metenp_basic_param2, 
    PPI_optstr, location_optstr, metlist_optstr].join("");

  //console.log(request_url);

  const req = await fetch(request_url);
  const res = await req.json();

  let MetENPobj = {RESTURL: request_url, res:res}
  return MetENPobj;
//  return res;
})
.story(props => ({
  abstract: `For the given list of metabolites, MetENP analysis is carried out and the result (a list of output files) is displayed as a table for further exploration.`,
  legend: `The list of results from the MetNet pipeline (combination of a subset of MetENP pipeline and an additional step of developing PPI network using STRING-DB APIson the list of genes corresponding to the metabolic enzymes). The components of the results can be accessed following the URLs for the corresponding text-based descriptions.`,
  introduction: `Given a list of metabolites, e.g., metabolites with significant change between two conditions such as disease/normal or treatment/control in a metabolomics study obtained by using MetENP\\ref{https://www.biorxiv.org/content/10.1101/2020.11.20.391912v1} or another tool, a researcher may want to find what are the pathways and functions affected. MetENP provides a list of metabolic pathways affected, but it does not provide such information for signaling/regulatory pathways affected, e.g., using protein-protein interaction (PPI) information. MetNet, which is a metabolite-centric tool, fills that gap. The user is strongly encouraged to understand MetENP before using MetNet. MetNet uses several features and code-base from MetENP such as metabolite name harmonization using REFMET, metabolite class enrichment, metabolic pathway enrichment and visualization and identification of reactions related to the given metabolites and the genes (coding for enzymes) catalyzing these reactions. The list of genes is used to develop their PPI network using STRING-DB APIs. Currently, MetNet is available as a REST API.`,
  methods: `Given a list of metabolites, MetNet uses several features and code-base from MetENP\\ref{https://www.biorxiv.org/content/10.1101/2020.11.20.391912v1} such as metabolite name harmonization using REFMET, metabolite class enrichment, metabolic pathway enrichment and visualization and identification of reactions related to the given metabolites and the genes (coding for enzymes) catalyzing these reactions. The list of genes is used to develop their PPI network using STRING-DB APIs. The results are provided as a json object or a table with the columns FileDescription and FileURL. The components of the results can be accessed following these URLs. Please note that in the REST interface, the order of the parameters is fixed.`,
}))
.build()

