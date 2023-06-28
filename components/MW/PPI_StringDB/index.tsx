import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { GeneTerm } from '@/components/core/input/term'
import { GeneSet } from '@/components/core/input/set'
import { GraphPlot } from '@/components/viz/graph'
import { metabolomicsworkbench_icon } from '@/icons'
import { z } from 'zod'
import { additional_info_icon, gene_icon } from '@/icons'
//import { VALID_LOADERS } from 'next/dist/shared/lib/image-config'


// How the schema validation works: https://codex.so/zod-validation-en

// StringDB PPI edge schema informaTION: FROM: https://string-db.org/cgi/help.pl?subpage=api%23getting-all-the-string-interaction-partners-of-the-protein-set
/*
Output fields (TSV and JSON formats):
Field 	Description
stringId_A 	STRING identifier (protein A)
stringId_B 	STRING identifier (protein B)
preferredName_A 	common protein name (protein A)
preferredName_B 	common protein name (protein B)
ncbiTaxonId 	NCBI taxon identifier
score 	combined score
nscore 	gene neighborhood score
fscore 	gene fusion score
pscore 	phylogenetic profile score
ascore 	coexpression score
escore 	experimental score
dscore 	database score
tscore 	textmining score
*/

export const StringDBedgeC = z.object({
    stringId_A: z.string().optional(),
    stringId_B: z.string().optional(),
    preferredName_A: z.string(),
    preferredName_B: z.string(),
    ncbiTaxonId: z.number(),
    score: z.number().optional(),
    nscore: z.number().optional(),
    fscore: z.number().optional(),
    pscore: z.number().optional(),
    ascore: z.number().optional(),
    escore: z.number().optional(),
    dscore: z.number().optional(),
    tscore: z.number().optional()
    })

    export const MyedgeC = z.object({
        SYMBOL_A: z.string(),
        SYMBOL_B: z.string(),
        score: z.number().optional()
    })

    export type StringDBedge = z.infer<typeof StringDBedgeC>
    export type Myedge = z.infer<typeof MyedgeC>

// important ref: https://rsinohara.github.io/json-to-zod-react/
export const StringDBedgeArrayC = z.array(
    StringDBedgeC
)
export const MyedgeArrayC = z.array(
    MyedgeC
)

export type StringDBedgeArray = z.infer<typeof StringDBedgeArrayC>
export type MyedgeArray = z.infer<typeof MyedgeArrayC>

export function SimplifyStringDBedgeArray(data:StringDBedgeArray) {
  // Given the list of edges in standard StringDB PPI format, extract only three columns: 
  // preferredName_A [call SYMBOL_A], preferredName_B [call SYMBOL_B] and score.
  let edges : MyedgeArray = data.map(d => {
   return {
    SYMBOL_A : d.preferredName_A, 
    SYMBOL_B : d.preferredName_B, 
    score : d.score,
  };
  });
  return edges;
}

export function SimplifyStringDBedge(d:StringDBedge) {
    // Given the edge in standard StringDB PPI format, extract only three columns: 
    // preferredName_A [call SYMBOL_A], preferredName_B [call SYMBOL_B] and score.
    let edge: Myedge =  { SYMBOL_A : d.preferredName_A, SYMBOL_B : d.preferredName_B, score : d.score};
    return edge;
}
  
export function GetAllNodes_from_StringDBedgeArray(data:StringDBedgeArray) {
    // Given the list of edges in standard StringDB PPI format, extract the array of all nodes (unique)
    let sources : string[] = Array.from(new Set( data.map(a => a.preferredName_A)));
    let targets : string[] = Array.from(new Set( data.map(a => a.preferredName_B)));
    let allnodes = [...new Set( sources.concat(targets) )] ;
    return allnodes;
}

export function GetAllNodes_from_MyedgeArray(data:MyedgeArray) {
    // Given the list of edges in MyedgeArray, extract the array of all nodes (unique)
    let sources : string[] = Array.from(new Set( data.map(a => a.SYMBOL_A)));
    let targets : string[] = Array.from(new Set( data.map(a => a.SYMBOL_B)));
    let allnodes = [...new Set( sources.concat(targets) )] ;
    return allnodes;
}

export function Format_StringDBedgeArray_for_GraphPlot(data:StringDBedgeArray) {
    // Given the list of edges in standard StringDB PPI format, get an object in the codec format for network plot:
    // https://github.com/nih-cfde/playbook-partnership/blob/network-viz/components/viz/graph/index.tsx
/*
.codec(z.object({
    nodes: z.array(z.object({
      id: z.string(),
      label: z.string().optional(),
      type: z.string(),
    })),
    edges: z.array(z.object({
      source: z.string(),
      target: z.string(),
    })),
  }))
*/
    let allnodes = GetAllNodes_from_StringDBedgeArray(data);
    let nodes_array = allnodes.map(d => {
        return {
         id : d, 
         label : d, 
         type : "gene",
       };
    });

    let edges_array = data.map(d => {
        return {
         source : d.preferredName_A, 
         target : d.preferredName_B, 
       };
    });

    return {nodes: nodes_array, edges: edges_array};
}
  
// A unique name for your data type is used here
export const StringDB_PPI_Network = MetaNode('StringDB_PPI_Network')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'StringDB PPI network',
    description: 'StringDB PPI network',
    icon: [gene_icon, additional_info_icon],
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using zod, a compile-time and runtime type-safe codec can be constructed
  .codec(StringDBedgeArrayC)
  // react component rendering your data goes here
  .view(data => {
    const dataobj = data; //const dataobj = [data][0];

    return(
          <div><table>
    <tr>
    <th>SYMBOL A</th>
    <th>SYMBOL B</th>
    <th>Score</th>
    </tr>
{dataobj.map((val, key) => {
      return (
      <tr key={key}>
        <td>{val.preferredName_A}</td>
        <td>{val.preferredName_B}</td>
        <td>{val.score}</td>
      </tr>
    )
    })}
    </table></div>
    )
  })
  .build()

  // A unique name for your resolver is used here
export const FetchStringDBPPI_Gene = MetaNode('FetchStringDBPPI')
// Human readble descriptors about this node should go here
.meta({
  label: 'Fetch StringDB PPI',
  description: 'Given a gene or gene set (SYMBOL), extract PPI using StringDB APIs',
  icon: [metabolomicsworkbench_icon],
})
// This should be a mapping from argument name to argument type
//  the types are previously defined Meta Node Data Types
.inputs({ gene: GeneTerm })
// This should be a single Meta Node Data Type
.output(StringDB_PPI_Network)
// The resolve function uses the inputs and returns output
//  both in the shape prescribed by the data type codecs
.resolve(async (props) => {
  // https://string-db.org/help/api/
  // Based on the API at: https://string-db.org/cgi/help.pl?subpage=api%23getting-all-the-string-interaction-partners-of-the-protein-set
  // See also: get, post, fetch, request: https://medium.com/meta-box/how-to-send-get-and-post-requests-with-javascript-fetch-api-d0685b7ee6ed
  const species_txid = "9606"
  const string_api_url = "https://version-11-5.string-db.org/api"
  const output_format = "json"
  const method = "interaction_partners"
  
  const request_url_base = [string_api_url, output_format, method].join("/");
  const my_genes = [props.inputs.gene];
  
  const params = {
      identifiers : my_genes.join("%0d"), // your protein
      species : species_txid, // species NCBI identifier 
      limit : "5000",
      required_score : "400",
      caller_identity : "sc-cfdewebdev.sdsc.edu" // your app name
  }
  const params_str = "identifiers=" + params.identifiers + "&" +
  "species=" + params.species + "&" +
  "limit=" + params.limit + "&" +
  "required_score=" + params.required_score + "&" +
  "caller_identity=" + params.caller_identity;

  const req = await fetch(`${request_url_base}?${params_str}`);
  const res = await req.json()

  return res
})
.story(props =>
  `For the given gene ID (SYMBOL), StringDB PPI was extracted using their API [\\ref{STRING api, https://string-db.org/cgi/help.pl?subpage=api%23getting-all-the-string-interaction-partners-of-the-protein-set}].`
)
.build()

// A unique name for your resolver is used here
export const StringDBPPI_to_GeneSet = MetaNode('StringDBPPI_to_GeneSet')
// Human readble descriptors about this node should go here
.meta({
  label: 'Given StringDB PPI, generate the list of nodes (GeneSet)',
  description: 'Given StringDB PPI, generate the list of nodes (GeneSet)',
  icon: [metabolomicsworkbench_icon],
})
// This should be a mapping from argument name to argument type
//  the types are previously defined Meta Node Data Types
.inputs({ data: StringDB_PPI_Network })
// This should be a single Meta Node Data Type
.output(GeneSet)
// The resolve function uses the inputs and returns output
//  both in the shape prescribed by the data type codecs
.resolve(async (props) => {
  const allnodes = GetAllNodes_from_StringDBedgeArray(props.inputs.data);  
  return {"description":"", "set":allnodes} ;
})
.story(props =>
  `For the Given StringDB PPI, the list of nodes (GeneSet) is generated.`
)
.build()

// A unique name for your resolver is used here
export const StringDBPPI_to_GraphPlot = MetaNode('StringDBPPI_to_GraphPlot')
// Human readble descriptors about this node should go here
.meta({
  label: 'Reformat StringDB PPI for plotting',
  description: 'Given StringDB PPI, reformat for plotting',
  icon: [metabolomicsworkbench_icon],
})
// This should be a mapping from argument name to argument type
//  the types are previously defined Meta Node Data Types
.inputs({ data: StringDB_PPI_Network })
// This should be a single Meta Node Data Type
.output(GraphPlot)
// The resolve function uses the inputs and returns output
//  both in the shape prescribed by the data type codecs
.resolve(async (props) => {
  const GraphPlotObj = Format_StringDBedgeArray_for_GraphPlot(props.inputs.data);

  return GraphPlotObj;
})
.story(props =>
  `For the Given StringDB PPI, the list of nodes (GeneSet) is generated.`
)
.build()
