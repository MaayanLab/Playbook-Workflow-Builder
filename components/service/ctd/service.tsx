import { MetaNode } from '@/spec/metanode'
import { FileURL } from '@/components/core/file'
import {CTDAdjacencyAndExpressions, CTD_DataSet} from './utils'
import { GeneSet } from '@/components/core/set'
import { z } from 'zod'
import { file_transfer_icon, datafile_icon, ctd_icon } from '@/icons'
import { fileAsStream } from  '@/components/core/file/api/download'
import { GraphPlot } from '@/components/viz/graph'
import { fileFromStream, uploadFile } from  '@/components/core/file/api/upload'
import FormData from 'form-data'
import { Readable } from 'stream'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { pythonStream } from '@/utils/python'


const CTDResponseC = z.object({
  "highlyConnectedGenes": z.any(),
  "guiltyByAssociationGenes": z.any(),
  "jsonGraph": z.object({
    "nodes": z.array(z.object({
      "id": z.string().optional(),
      "name": z.string().optional(),
      "type": z.string(),
      "color": z.string().optional()
    })),
    "edges": z.any().optional(),
    "interactions": z.array(z.object({
      "source": z.string(),
      "target": z.string()
    })).optional(),
  }).nullable(),
  "report": z.object({
    "type": z.string(),
    "message": z.any().optional()
  }).nullable()
});
export type CTDResponse = z.infer<typeof CTDResponseC>

export async function getCTDGenSetResponse(strValue: string): Promise<CTDResponse> {
  const res = await fetch(`http://genboree.org/pb-ctd/rest/playbook_ctd/ctd/geneList`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
    },
    body: strValue
  })
  return await res.json()
}

export async function getCTDFileResponse(formData: FormData): Promise<CTDResponse> {
  const { default: axios } = await import('axios')
  const res = await axios.post(`http://genboree.org/pb-ctd/rest/playbook_ctd/ctd/file`, formData, {
    headers: { ...formData.getHeaders() },
    responseType: 'json',
  })
  return res.data
}

export async function getCTDAdjacency(formData: FormData): Promise<Readable> {
  const { default: axios } = await import('axios')
  const res = await axios.post(`http://genboree.org/pb-ctd/rest/playbook_ctd/ctd/createCustomMatrix`, formData, {
    headers: { ...formData.getHeaders() },
    responseType: 'stream'
  });
  if (res.status !== 200) throw new Error(await res.data.read())
  return res.data
}

export async function getCTDPrecalculations(formData: FormData): Promise<Readable> {
  const { default: axios } = await import('axios')
  const res = await axios.post(`http://genboree.org/pb-ctd/rest/playbook_ctd/ctd/getCustomPermutations`, formData, {
    headers: { ...formData.getHeaders() },
    responseType: 'stream',
  });
  if (res.status !== 200) throw new Error(await res.data.read())
  return res.data
}

export async function getCTDCustomResponse(formData: FormData): Promise<CTDResponse> {
  const { default: axios } = await import('axios')
  const res = await axios.post(`http://genboree.org/pb-ctd/rest/playbook_ctd/ctd/useCustomMatrix`, formData, {
    headers: { ...formData.getHeaders() },
    responseType: 'json'
  });
  return res.data
}

  export const CTD_CreateACustomMatrix = MetaNode('CTD_CreateACustomMatrix')
  .meta({
    label: `CTD - Create a Custom Adjacency File`,
    description: "Create a Custom Adjacency JSON File Using a Gene Expressions File.",
    tags: {
      Type: {
        CTD: 1
      }
    },
    icon: [ctd_icon],
    external: true,
  })
  .inputs({geneExpressions: GeneCountMatrix})
  .output(CTDAdjacencyAndExpressions)
  .resolve(async (props) => {
    const fileReader = pythonStream('components.service.ctd.csv_read_stream', {
      kargs: [props.inputs.geneExpressions],
    })

    const formData = new FormData();
    formData.append('csvExpressionsFile', fileReader, props.inputs.geneExpressions.filename);

    const response = await getCTDAdjacency(formData);
    if(response == null){
      throw new Error("Adjacency file is null or empty!");
    }
    const adjacencyFile = await uploadFile(await fileFromStream(response, 'ctdCustomAdjacency.json'));
    
    let output = {
      'adjacency':adjacencyFile,
      'expressions': props.inputs.geneExpressions
    }
    return output;
  }).story(props => ({
    abstract: `A custom CTD Adjacency JSON File is being created!`,
    introduction: `Connect the Dots (CTD) is an algorithm that was developed in order to quickly identify highly connected subsets of graph nodes\\ref{doi:10.1371/journal.pcbi.1009551}.`,
    methods: `A custom CTD Adjacency JSON File is created for the input gene count matrix file by using a CTD adjacency JSON file and gene expression CSV file\\ref{doi:10.1371/journal.pcbi.1009551}.`,
    tableLegend: `A table displaying highly connected sets of proteins from the input graph.`,
  })).build()

export const Execute_CTD_Precalculations_Combined = MetaNode('Execute_CTD_Precalculations_Combined')
  .meta({
    label: `Connect the Dots in Precalculated Graph`,
    description: `Use CTD to "Connect the Dots" and identify highly connected set of proteins using the pre-calculated graph.`,
    icon: [ctd_icon],
    external: true,
  })
  .inputs({ geneSet: GeneSet, ctdAdjacencyAndExpressions: CTDAdjacencyAndExpressions})
  .output(CTD_DataSet)
  .resolve(async (props) => {
    let geneNamesList = props.inputs.geneSet.set;
    let ctdAdjacencyFile = props.inputs.ctdAdjacencyAndExpressions.adjacency;
    const ctdAdjacencyFileReader = await fileAsStream(ctdAdjacencyFile);


    console.log("Adjacency file: "+ctdAdjacencyFile.filename);
    console.log("Expressions file: "+props.inputs.ctdAdjacencyAndExpressions.expressions.filename);

    const formData = new FormData(); 
    formData.append('geneList', geneNamesList.join('\n'), { filename: 'geneSetTempFile.csv', contentType: 'text/plain' })
    formData.append('customAdjacency', ctdAdjacencyFileReader, ctdAdjacencyFile.filename);

    const response = await getCTDPrecalculations(formData);
    if(response == null){
      throw new Error("Permutations file is null or empty!");
    }
    const permutationsfile = await uploadFile(await fileFromStream(response, "ctdPermutations.json"));
    
    let output = {
      'ctdPermutations':permutationsfile,
      'expressions': props.inputs.ctdAdjacencyAndExpressions.expressions,
      'geneSet': geneNamesList
    }
    return output;
  }).story(props => ({
    abstract: `Input Gene Set and Adj. Matrix to send to the CTD API for precalculations.`,
    introduction: `Connect the Dots (CTD) is an algorithm that was developed in order to quickly identify highly connected subsets of graph nodes\\ref{doi:10.1371/journal.pcbi.1009551}.`,
    methods: `Highly connected sets of proteins are identified using the CTD algorithm and a precalculated graph including a gene set and adjacency matrix\\ref{doi:10.1371/journal.pcbi.1009551}.`,
    tableLegend: `A table displaying highly connected sets of proteins from the input graph.`,
  })).build()

  export const CTDResponseInfo = MetaNode('CTDResponseInfo')
  .meta({
    label: 'CTD Response Info',
    description: '',
    icon: [datafile_icon]
  })
  .codec(CTDResponseC)
  .view(props => {

    let message = "";
    if(props.report != null && props.report.type == 'note'){
      message = "Note: "+props.report.message;
    }

    let hcgNum = 0;
    if(props.highlyConnectedGenes != null){
      hcgNum = props.highlyConnectedGenes.length;
    }

    let gbaNum = 0;
    if(props.guiltyByAssociationGenes != null){
      gbaNum = props.guiltyByAssociationGenes.length;
    }

    let nodesNum = 0;
    if(props.jsonGraph != null){
      nodesNum = props.jsonGraph.nodes.length;
    }

    return(
      <div className="prose max-w-none">
        <p>{message}</p>
        <p>Highly Connected Genes num: {hcgNum}</p>
        <p>Guilty By Association Genes num: {gbaNum}</p>
        <p>Json Graph Nodes num: {nodesNum}</p>
      </div>
    )
  }).build()

  export const CTD_UseCustomMatrixCombined = MetaNode('CTD_UseCustomMatrixCombined')
  .meta({
    label: `CTD Custom Response - Final`,
    description: `This card sends the gene expression matrix, the gene set and the adjacency matrix to the CTD API for processing. The following card shows the outputs for CTD.`,
    external: true,
  })
  .inputs({ ctdDataSet: CTD_DataSet})
  .output(CTDResponseInfo)
  .resolve(async (props) => {
    let geneNamesList = props.inputs.ctdDataSet.geneSet;
    let expressionsFile = props.inputs.ctdDataSet.expressions;
    let permutationsFile = props.inputs.ctdDataSet.ctdPermutations;

    const expressionsFileReader = await fileAsStream(expressionsFile);
    const permutationsFileReader = await fileAsStream(permutationsFile);
    const formData = new FormData();
    formData.append('csvGenesFile', geneNamesList.join('\n'), { filename: 'geneSetTempFile.csv', contentType: 'text/plain' });
    formData.append('csvExpressionsFile', expressionsFileReader, expressionsFile.filename);
    formData.append('jsonPermutationsFile', permutationsFileReader, permutationsFile.filename);

    let response = await getCTDCustomResponse(formData);
    if(response != null && response.report != null && response.report.type == 'error'){
      throw new Error(response.report.message);
    }
   
    return response;
  }).story(props => ({
    abstract: `The three files where send to the CTD API for precalculations.`,
    introduction: `Connect the Dots (CTD) is an algorithm that was developed in order to quickly identify highly connected subsets of graph nodes\\ref{doi:10.1371/journal.pcbi.1009551}.`,
    methods: `The custom gene list, adjacency matrix file and permutations (RData) file were sent to the CTD API for precalculations.`,
    tableLegend: `A table displaying the CTD output including highly connected node subsets.`,
  })).build()

  export const Highly_Connected_Genes = MetaNode('Highly_Connected_Genes')
  .meta({
    label: `Extract Highly Connected Genes`,
    description: "Extract nodes that are determined to be highly connected by CTD in your initial node set and display them in a table."
  })
  .inputs({ ctdResponseInfo: CTDResponseInfo })
  .output(GeneSet)
  .resolve(async (props) => {
    let highlyConnectedGenes = props.inputs.ctdResponseInfo.highlyConnectedGenes;

    if(highlyConnectedGenes == null || highlyConnectedGenes.length == 0){
      throw new Error("No Highly Connected Genes available, please use a different input gene set!");
    }

    let geneSet = {
      description: 'Highly_Connected_Genes, CTD response',
      set:  highlyConnectedGenes
    };
    return geneSet;
  }).story(props => ({
    abstract: `A list of Highly Connected Genes was obtained from the CTD output.`,
    introduction: `Connect the Dots (CTD) is an algorithm that was developed in order to quickly identify highly connected subsets of graph nodes\\ref{doi:10.1371/journal.pcbi.1009551}.`,
    methods: `Nodes within the initial node set are extracted if they are determined to be highly connected by CTD. `,
    tableLegend: `A table of the extracted highly connected genes as determined by CTD.`,
  })).build()

export const Guilty_By_Association_Genes = MetaNode('Guilty_By_Association_Genes')
  .meta({
    label: `Extract Guilty By Association Genes`,
    description: `Extract nodes that are "guilty by association" and connect your initial genes of interest within the graph.`
  })
  .inputs({ ctdResponseInfo: CTDResponseInfo })
  .output(GeneSet)
  .resolve(async (props) => {
    let guiltyByAssociationGenes = props.inputs.ctdResponseInfo.guiltyByAssociationGenes;

    if(guiltyByAssociationGenes == null || guiltyByAssociationGenes.length == 0){
      throw new Error("No Guilty By Association Genes available, please use a different input gene set!");
    }

    let geneSet = {
      description: 'Guilty_By_Association_Genes, CTD response',
      set:  guiltyByAssociationGenes
    };
    return geneSet;
  }).story(props => ({
    abstract: `A list of Guilty By Association Genes was obtained from the CTD output.`,
    introduction: `Connect the Dots (CTD) is an algorithm that was developed in order to quickly identify highly connected subsets of graph nodes\\ref{doi:10.1371/journal.pcbi.1009551}.`,
    methods: `Nodes are extracted as Guilty By Association if they connet initial genes of interest within the graph.`,
    tableLegend: `A table of Guilty By Association genes.`,
  })).build()

export const CTD_Graph_Nodes = MetaNode('CTD_Graph_Nodes')
  .meta({
    label: 'CTD Graph',
    description: `This is a visual display of the nodes that are found to be highly connected with CTD. Nodes that connect these nodes and are "guilty by association" will also be displayed.`,
  })
  .inputs({ ctdResponseInfo: CTDResponseInfo })
  .output(GraphPlot)
  .resolve(async (props) => {
    let jsonGraph = null;
    if(props.inputs.ctdResponseInfo.jsonGraph != null){
      jsonGraph = props.inputs.ctdResponseInfo.jsonGraph;
    }

    if(jsonGraph == null || jsonGraph.nodes == null || jsonGraph.nodes.length == 0){
      throw new Error("No Gene Graph Nodes available, please use a different input gene set!");
    }

    let highlyConnectedGenes = props.inputs.ctdResponseInfo.highlyConnectedGenes;
    let hcN = highlyConnectedGenes.length;
    let guiltyByAssociationGenes = props.inputs.ctdResponseInfo.guiltyByAssociationGenes;
    let gbaN = guiltyByAssociationGenes.length;
    let n = 0;
    if(hcN > gbaN){
      n = hcN;
    }else{
      n = gbaN
    }

    let graphNodes = jsonGraph.nodes;  
    for(let i=0; i<n; i++){
      let hcGene = null;
      try {
        hcGene = highlyConnectedGenes[i];
      } catch (error) {}      
      
      let gbaGene = null;
      try {
        gbaGene = guiltyByAssociationGenes[i];  
      } catch (error) {}

      let chChecked = false;
      let gbaChecked = false;
      for(let k in graphNodes){
        if(!chChecked && hcGene != null && hcGene == graphNodes[k].name){
          graphNodes[k].type = "hc";
          graphNodes[k].color = "#ff8566";
          chChecked = true;
        }else if(!gbaChecked && gbaGene != null && gbaGene == graphNodes[k].name){
          graphNodes[k].type = "gba";
          graphNodes[k].color = "#99bbff";
          gbaChecked = true;
        }

        if(chChecked && gbaChecked){
          break;
        }
      }

    }

    return {
      edges: jsonGraph.interactions || [],
      nodes: graphNodes.map(({ name, ...rest }) => ({ ...rest, id: name || '(no name)' })),
    }
  }).story(props => ({
    abstract: `Graph Nodes were extracted from the CTD output.`,
    introduction: `Connect the Dots (CTD) is an algorithm that was developed in order to quickly identify highly connected subsets of graph nodes\\ref{doi:10.1371/journal.pcbi.1009551}.`,
    methods: `Highly connected nodes as well as nodes that are connect these nodes (Guilty By Association) are extracted.`,
    tableLegend: `A table of graph nodes extracted using CTD.`,
  })).build()

  export const GeneSet_CTD_String = MetaNode('GeneSet_CTD_String')
  .meta({
    label: `CTD - Connect the Dots in STRING`,
    description: 'Use CTD to "Connect the Dots" and identify highly connected set of proteins in the STRING protein interaction graph. *Please note 10-150 genes of interest are required to run CTD',
    icon: [ctd_icon],
    external: true,
  })
  .inputs({ geneset: GeneSet })
  .output(CTDResponseInfo)
  .resolve(async (props) => {
    console.log("CTD String, gen set processing.");
    let requestBody = {
      "graphType": "string",
      "geneList": props.inputs.geneset.set
    }

    let response = await getCTDGenSetResponse(JSON.stringify(requestBody));
    if(response != null && response.report != null && response.report.type == 'error'){
      throw new Error(response.report.message);
    }

    return response;
  }).story(props => ({
    abstract: `CTD is applied which diffuses through all nodes in STRING\\ref{doi:10.1093/nar/gku1003} to identify nodes that are "guilty by association" and highly connected to the initial gene set of interest\\ref{doi:10.1371/journal.pcbi.1009551}\\ref{doi:10.1016/j.isci.2022.105799}.`,
    introduction: `Connect the Dots (CTD) is an algorithm that was developed in order to quickly identify highly connected subsets of graph nodes\\ref{doi:10.1371/journal.pcbi.1009551}. STRING is an online database that stores graphical protein-protein interaction networks\\ref{doi:10.1093/nar/gku1003}.`,
    methods: `CTD is used to identify the highly connected and Guilty By Association sets of proteins in the STRING protein interaction graph\\ref{doi:10.1093/nar/gku1003}.`,
    tableLegend: `A table displaying the CTD output including highly connected node subsets.`,
  })).build()

export const GeneSet_CTD_Wikipathways = MetaNode('GeneSet_CTD_Wikipathways')
  .meta({
    label: `CTD - Connect the Dots in Wikipathways`,
    description: `Use CTD to "Connect the Dots" and identify highly connected set of genes in the WikiPathways pathway annotation graph. *Please note 10-150 genes of interest are required to run CTD`,
    icon: [ctd_icon],
    external: true,
  })
  .inputs({ geneset: GeneSet })
  .output(CTDResponseInfo)
  .resolve(async (props) => {
    console.log("CTD Wikipathways, gen set processing.");
    let requestBody = {
      "graphType": "wikipathways",
      "geneList": props.inputs.geneset.set
    }

    let response = await getCTDGenSetResponse(JSON.stringify(requestBody));
    if(response != null && response.report != null && response.report.type == 'error'){
      throw new Error(response.report.message);
    }
    return response;
  }).story(props => ({
    abstract: `CTD is applied which diffuses through all nodes in WikiPathways\\ref{doi:10.1093/nar/gkad960} to identify nodes that are "guilty by association" and highly connected to the initial gene set of interest\\ref{doi:10.1371/journal.pcbi.1009551}\\ref{doi:10.1016/j.isci.2022.105799}.`,
    introduction: `Connect the Dots (CTD) is an algorithm that was developed in order to quickly identify highly connected subsets of graph nodes\\ref{doi:10.1371/journal.pcbi.1009551}. WikiPathways is a web-based software application that stores biological pathways contributed by the research community\\ref{doi:10.1093/nar/gkad960}. `,
    methods: `CTD is used to identify the highly connected and Guilty By Association nodes in WikiPathways\\ref{doi:10.1093/nar/gkad960}\\ref{doi:10.1371/journal.pcbi.1009551}.`,
     tableLegend: `A table displaying the CTD output including highly connected node subsets.`,
   })).build()

export const GenesFile_CTD_String = MetaNode('GenesFile_CTD_String')
  .meta({
    label: `CTD String For Genes Set File`,
    description: "Ensure a file contains a gene set, values separated by a \\n character  and with the extension .csv",
    icon: [file_transfer_icon],
    external: true,
  })
  .inputs({ file: FileURL })
  .output(CTDResponseInfo)
  .resolve(async (props) => {
    const fileReader = await fileAsStream(props.inputs.file);
    console.log("CTD String, file processing: "+props.inputs.file.filename);
    const formData = new FormData();
    formData.append('csvGenesFile', fileReader, props.inputs.file.filename);
    formData.append('graphType', "string");

    let response = await getCTDFileResponse(formData);
    if(response != null && response.report != null && response.report.type == 'error'){
      throw new Error(response.report.message);
    }
    return response;
  })
  .story(props => ({
    abstract: `CTD is applied which diffuses through all nodes in WikiPathways\\ref{doi:10.1093/nar/gkad960}] to identify nodes that are "guilty by association" and highly connected to the initial gene set of interest\\ref{doi:10.1371/journal.pcbi.1009551}\\ref{doi:10.1016/j.isci.2022.105799}.`,
    introduction: `Connect the Dots (CTD) is an algorithm that was developed in order to quickly identify highly connected subsets of graph nodes\ref{doi:10.1371/journal.pcbi.1009551}.`,
    methods: `Input files were checked for proper formatting before being used by CTD. Proper file formatting is defined as having values separated by a \\n character and a .csv extension.`,
    tableLegend: `A table displaying the CTD output including highly connected node subsets.`,
  }))
  .build()

  export const GenesFile_CTD_Wikipathways = MetaNode('GenesFile_CTD_Wikipathways')
  .meta({
    label: `CTD Wikipathways For Genes Set File`,
    description: `Ensure a file contains a gene set, values separated by a \\n character and with the extension .csv`,
    icon: [file_transfer_icon],
    external: true,
  })
  .inputs({ file: FileURL })
  .output(CTDResponseInfo)
  .resolve(async (props) => {
    const fileReader = await fileAsStream(props.inputs.file);
    console.log("CTD Wikpathways, file processing: "+props.inputs.file.filename);
    const formData = new FormData();
    formData.append('csvGenesFile', fileReader, props.inputs.file.filename);
    formData.append('graphType', "wikipathways");

    let response = await getCTDFileResponse(formData);
    if(response != null && response.report != null && response.report.type == 'error'){
      throw new Error(response.report.message);
    }
    return response;
  })
  .story(props => ({
    abstract: `CTD is applied which diffuses through all nodes in WikiPathways\\ref{doi:10.1093/nar/gkad960}] to identify nodes that are "guilty by association" and highly connected to the initial gene set of interest\\ref{doi:10.1371/journal.pcbi.1009551}\\ref{doi:10.1016/j.isci.2022.105799}.`,
    introduction: `Connect the Dots (CTD) is an algorithm that was developed in order to quickly identify highly connected subsets of graph nodes\ref{doi:10.1371/journal.pcbi.1009551}.`,
    methods: `Input files were checked for proper formatting before being used by CTD. Proper file formatting is defined as having values separated by a \\n character and a .csv extension.`,
    tableLegend: `A table displaying the CTD output including highly connected node subsets.`,
  })).build()
