import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { GeneTerm } from '@/components/core/term'
import { GeneSet } from '@/components/core/set'
import { GraphPlot } from '@/components/viz/graph'
import { metabolomicsworkbench_icon } from '@/icons'
import { z } from 'zod'
import { additional_info_icon, gene_icon } from '@/icons'
import * as array from '@/utils/array'
import { GetGeneSetIDConv } from '@/components/MW/ConvertedGeneID'

// How the schema validation works: https://codex.so/zod-validation-en

// StringDB PPI edge schema information: FROM: https://string-db.org/cgi/help.pl?subpage=api%23getting-all-the-string-interaction-partners-of-the-protein-set
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

export const StringDB_PPI_NetworkDataC = z.object({
  species_txid: z.string(),
  species_id: z.string(),
  edges: StringDBedgeArrayC
})

export type StringDB_PPI_NetworkData = z.infer<typeof StringDB_PPI_NetworkDataC>

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
    return array.unique(data.flatMap(a => [a.preferredName_A, a.preferredName_B]));
}

export function GetAllNodes_from_MyedgeArray(data:MyedgeArray) {
    // Given the list of edges in MyedgeArray, extract the array of all nodes (unique)
    return array.unique(data.flatMap(a => [a.SYMBOL_A, a.SYMBOL_B]));
}

export async function Format_StringDBedgeArray_for_GraphPlot(data:StringDBedgeArray, species_id:string = "hsa", geneid_type:string = "SYMBOL_OR_ALIAS") {
    // Need to use async since GetGeneSetIDConv is an async function
    // Given the list of edges in standard StringDB PPI format, get an object in the codec format for network plot:
    // https://github.com/MaayanLab/Playbook-Workflow-Builder/blob/network-viz/components/viz/graph/index.tsx
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
    let allnodes_more = await GetGeneSetIDConv(allnodes,  species_id, geneid_type); // GetGeneSetIDConv is async function, so, await is needed
    let nodes_array = allnodes_more.map((d,i) => {
        let ezid = d.ENTREZID;
        return {
         id :  d.SYMBOL, 
         label : d.SYMBOL, // allnodes_more[i].GENENAME, // d, // need to cross-check order if d goes with allnodes_more[i]
         link: `https://www.ncbi.nlm.nih.gov/gene/${ezid}`,
         hovertext: d.GENENAME,
         type : (d.SYMBOL==allnodes[0]) ? "gene1" : "gene",
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
  .codec(StringDB_PPI_NetworkDataC)
  // react component rendering your data goes here
  .view(data => {
    const dataobj = data.edges; //const dataobj = [data][0];

    return(
          <div className="prose max-w-none">
            <h3>Species: {data.species_id}</h3>
            <h3>Taxonomy ID: {data.species_txid}</h3>
            <table>
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
export const FetchStringDBPPI = MetaNode('FetchStringDBPPI')
// Human readble descriptors about this node should go here
.meta({
  label: 'Fetch StringDB PPI',
  description: 'Given a gene or gene set (SYMBOL), extract PPI using StringDB APIs',
  icon: [metabolomicsworkbench_icon],
  external: true,
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
  const species_id = "hsa"
  const string_api_url = "https://version-11-5.string-db.org/api"
  const output_format = "json"
  const method = "interaction_partners"
  
  const request_url_base = [string_api_url, output_format, method].join("/");
  const my_genes = [props.inputs.gene];
  
  // Try with the gene RPE as its PPI with 900 cut off is not that large and all those genes
  // don't have too many metabolites associated with them.
  const params = {
      identifiers : my_genes.join("%0d"), // your protein
      species : species_txid, // species NCBI identifier 
      limit : "5000",
      required_score : "900",
      caller_identity : "sc-cfdewebdev.sdsc.edu" // your app name
  }
  const params_str = "identifiers=" + params.identifiers + "&" +
  "species=" + params.species + "&" +
  "limit=" + params.limit + "&" +
  "required_score=" + params.required_score + "&" +
  "caller_identity=" + params.caller_identity;

  const req = await fetch(`${request_url_base}?${params_str}`);
  const res = await req.json()
  //return res
  let PPIobj = {species_id: species_id, species_txid: species_txid, edges:res};

  return PPIobj
})
.story(props => ({
  abstract: `For the given gene ID (SYMBOL), StringDB PPI was extracted using their API\\ref{doi:10.1093/nar/gkac1000}.`,
  introduction: `StringDB is a database of protein-protein interactions. Given a list of genes in terms of Gene Symbols, their API provides an easy way to fetch the list of protein-protein interactions\\ref{doi:10.1093/nar/gkac1000}.`,
  methods: `For ${props.inputs?.gene ? props.inputs.gene : 'the gene'}, StringDB\\ref{doi:10.1093/nar/gkac1000} PPI was extracted using their API.`,
  legend: `A StringDB\\ref{doi:10.1093/nar/gkac1000} PPI Network for ${props.inputs?.gene ? props.inputs.gene : 'the gene'}.`,
}))
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
  const allnodes = GetAllNodes_from_StringDBedgeArray(props.inputs.data.edges);  
  return {"description": "", "set": allnodes} ;
})
.story(props => ({
  abstract: `For the Given StringDB PPI, the list of nodes (GeneSet) is generated.`,
  introduction: `StringDB is a database of protein-protein interactions. Given a list of genes in terms of Gene Symbols, their API provides an easy way to fetch the list of protein-protein interactions\\ref{doi:10.1093/nar/gkac1000}.`,
  methods: `For the given StringDB\\ref{doi:10.1093/nar/gkac1000} PPI ${props.input_refs?.data}, the gene set is produced with the union of all nodes in the PPI network.`,
  legend: `A gene set corresponding to all unique genes in the StringDB\\ref{doi:10.1093/nar/gkac1000} PPI network from ${props.input_refs?.data}.`,
}))
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
  let species_id = props.inputs.data.species_id;
  let species_txid = props.inputs.data.species_txid;
  let edges = props.inputs.data.edges;
  
  const GraphPlotObj = Format_StringDBedgeArray_for_GraphPlot(edges, species_id);

  return GraphPlotObj;
})
.story(props => ({
  abstract: `For the Given StringDB PPI, the list of nodes (Gene Set) is generated.`,
  introduction: `StringDB is a database of protein-protein interactions. Given a list of genes in terms of Gene Symbols, their API provides an easy way to fetch the list of protein-protein interactions.`,
  methods: `For the Given StringDB PPI, an object is generated to be used with plotting the PPI network.`,
  legend: `PPI network for the list of genes submitted.`,
}))
.build()
