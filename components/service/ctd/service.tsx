import { MetaNode } from '@/spec/metanode'
import { FileURL, FileInput, FilePrompt, FileC } from '@/components/core/file'
import { CTDPrecalculationsFileURLs, CTDUseCustomMatrixFileURLs } from './input'
import { GeneSet } from '@/components/core/input/set'
import { z } from 'zod'
import { file_transfer_icon, datafile_icon } from '@/icons'
import { fileAsStream } from  '@/components/core/file/api/download'
import { GraphPlot } from '@/components/viz/graph'
import { fileFromStream, uploadFile } from  '@/components/core/file/api/upload'
import { downloadUrl } from '@/utils/download'
import { Table, Cell, Column} from '@/app/components/Table'
import FormData from 'form-data'
import { Readable } from 'stream'
import SafeRender from '@/utils/saferender'

const CTDResponseC = z.object({
  "highlyConnectedGenes": z.any(),
  "guiltyByAssociationGenes": z.any(),
  "jsonGraph": z.object({
    "nodes": z.array(z.object({
      "id": z.string().optional(),
      "name": z.string().optional(),
      "type": z.string()
    })),
    "edges": z.any().optional(),
    "interactions": z.array(z.object({
      "source": z.string(),
      "target": z.string()
    })).optional(),
  }),
  "report": z.object({
    "type": z.string(),
    "message": z.any().optional()
  })
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

export async function getCTDPrecalculationsResponse(formData: FormData): Promise<Readable> {
  const { default: axios } = await import('axios')
  const res = await axios.post(`http://genboree.org/pb-ctd/rest/playbook_ctd/getCtdCustomMatrix`, formData, {
    headers: { ...formData.getHeaders() },
    responseType: 'stream',
  });
  if (res.status !== 200) throw new Error(await res.data.read())
  return res.data
}

export async function getCTDUseCustomMatrix(formData: FormData): Promise<CTDResponse> {
  const { default: axios } = await import('axios')
  const res = await axios.post(`http://genboree.org/pb-ctd/rest/playbook_ctd/ctd/useCustomMatrix`, formData, {
    headers: { ...formData.getHeaders() },
    responseType: 'json'
  });
  return res.data
}

/*
export async function getCustomMatrixFromExpressions(formData: FormData): Promise<Readable> {
  const res = await axios.post(`http://genboree.org/pb-ctd/rest/playbook_ctd/ctd/createCustomMatrix`, formData, {
    headers: { ...formData.getHeaders() }
  });
  return res.data
}*/

export const CTD_FileDownload = MetaNode('CTD_FileDownload')
  .meta({
    label: 'CTD File Download',
    description: '',
    icon: [datafile_icon]
  })
  .codec(
    z.object({
      url: z.string(),
      filename: z.string(),
    })
  )
  .view(props => {
    return (
      <Table
      cellRendererDependencies={[1]}
      numRows={1}
      downloads={{
        'URL': () => downloadUrl(props.url, props.filename)
      }}
      >
      <Column
        name="Filename"
        cellRenderer={row => <Cell key={row+''}>{props.filename}</Cell>}
      />
      <Column
        name="URL"
        cellRenderer={row => <Cell key={row+''}>{props.url}</Cell>}
      />
      </Table>
    )
  }).build()


export const Execute_CTD_Precalculations = MetaNode('Execute_CTD_Precalculations')
  .meta({
    label: `CTD Precalculations With Custom Matrix`,
    description: 'Use CTD to "Connect the Dots" and identify highly connected set of proteins using the pre-calculated graph.'
  })
  .inputs({ ctdPrecalculationsFileURLs: CTDPrecalculationsFileURLs })
  .output(CTD_FileDownload)
  .resolve(async (props) => {
    let { geneListFile, adjMatrixFile } = props.inputs.ctdPrecalculationsFileURLs;

    const geneListFileReader = await fileAsStream(geneListFile);
    const adjMatrixFileReader = await fileAsStream(adjMatrixFile);
    const formData = new FormData();
    formData.append('geneList', geneListFileReader, geneListFile.filename);
    formData.append('customMatrix', adjMatrixFileReader, adjMatrixFile.filename);

    const response = await getCTDPrecalculationsResponse(formData);
    const file = await uploadFile(await fileFromStream(response, `derived.${"customRDataFile.RData"}`))
    return file;
  }).story(props =>   
    "The two files where send to the CTD API for precalculations."
  ).build()

export const Execute_CTD_Precalculations_Hybrid = MetaNode('Execute_CTD_Precalculations_Hybrid')
  .meta({
    label: `Connect the Dots in Precalculated Graph`,
    description: 'Use CTD to "Connect the Dots" and identify highly connected set of proteins using the pre-calculated graph.',
    hidden: true,
  })
  .inputs({ geneSet: GeneSet, file: FileURL})
  .output(CTD_FileDownload)
  .resolve(async (props) => {
    let geneNamesList = props.inputs.geneSet.set;
    let adjMatrixFile = props.inputs.file;
    const adjMatrixFileReader = await fileAsStream(adjMatrixFile);

    console.log("geneNamesList: "+JSON.stringify(geneNamesList))
    console.log("file: "+adjMatrixFile.filename)

    const formData = new FormData();
    formData.append('geneList', geneNamesList.join('\n'), { filename: 'geneSetTempFile.csv', contentType: 'text/plain' })
    formData.append('customMatrix', adjMatrixFileReader, adjMatrixFile.filename);

    const response = await getCTDPrecalculationsResponse(formData);
    const file = await uploadFile(await fileFromStream(response, `derived.${"customRDataFile.RData"}`))
    return file;
  }).story(props =>   
    "Input Gene Set and Adj. Matrix to send to the CTD API for precalculations."
  ).build()

export const Execute_CTD_Precalculations_Hybrid_Upload = MetaNode('Execute_CTD_Precalculations_Hybrid_Upload')
  .meta({
    label: `Connect the Dots in Precalculated Graph`,
    description: 'Use CTD to "Connect the Dots" and identify highly connected set of proteins using the pre-calculated graph.'
  })
  .codec(FileC)
  .inputs({ geneSet: GeneSet })
  .output(CTD_FileDownload)
  .prompt(props => <><FilePrompt {...props} />{props.output ? <SafeRender component={CTD_FileDownload.view} props={props.output} /> : null}</>)
  .resolve(async (props) => {
    let geneNamesList = props.inputs.geneSet.set;
    let adjMatrixFile = props.data;
    const adjMatrixFileReader = await fileAsStream(adjMatrixFile);

    console.log("geneNamesList: "+JSON.stringify(geneNamesList))
    console.log("file: "+adjMatrixFile.filename)

    const formData = new FormData();
    formData.append('geneList', geneNamesList.join('\n'), { filename: 'geneSetTempFile.csv', contentType: 'text/plain' })
    formData.append('customMatrix', adjMatrixFileReader, adjMatrixFile.filename);

    const response = await getCTDPrecalculationsResponse(formData);
    const file = await uploadFile(await fileFromStream(response, `derived.${"customRDataFile.RData"}`))
    return file;
  }).story(props =>   
    "Input Gene Set and Adj. Matrix to send to the CTD API for precalculations."
  ).build()

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

  /*
  export const CTD_CreateACustomMatrix = MetaNode('CTD_CreateACustomMatrix')
  .meta({
    label: `CTD Create a Custom Matrix`,
    description: "Create a Custom Matrix Using a Gene Expressions File."
  })
  .inputs({file: FileURL})
  .output(CTD_FileDownload)
  .resolve(async (props) => {
    const fileReader = await fileAsStream(props.inputs.file);
    console.log("CTD-Matrix from gene expressions file: "+props.inputs.file.filename);

    const formData = new FormData();
    formData.append('csvExpressionsFile', fileReader, props.inputs.file.filename);

    const response = await getCustomMatrixFromExpressions(formData);
    const file = await uploadFile(await fileFromStream(response, `derived.${"customMatrix.csv"}`))
    return file;
  }).story(props =>   
    "The three files where send to the CTD API for precalculations."
  ).build()
  */

  export const CTD_UseCustomMatrix = MetaNode('CTD_UseCustomMatrix')
  .meta({
    label: `CTD Response With Custom Matrix`,
    description: "Get CTD Reponse using a custom gene list, ajd. matrix file and RData file. Use the \"CTD Precalculations With Custom Matrix\" card to create the custom  RData file.!"
  })
  .inputs({ ctdUseCustomMatrixFileURLs: CTDUseCustomMatrixFileURLs })
  .output(CTDResponseInfo)
  .resolve(async (props) => {
    let { geneListFile, adjMatrixFile, rDataFile } = props.inputs.ctdUseCustomMatrixFileURLs;

    const geneListFileReader = await fileAsStream(geneListFile);
    const adjMatrixFileReader = await fileAsStream(adjMatrixFile);
    const rDataFileReader = await fileAsStream(rDataFile);
    const formData = new FormData();
    formData.append('csvGenesFile', geneListFileReader, geneListFile.filename);
    formData.append('customMatrix', adjMatrixFileReader, adjMatrixFile.filename);
    formData.append('customRData', rDataFileReader, rDataFile.filename);

    let response = await getCTDUseCustomMatrix(formData);
    if(response != null && response.report != null && response.report.type == 'error'){
      throw new Error(response.report.message);
    }
    return response;
  }).story(props =>   
    "The three files where send to the CTD API for precalculations."
  ).build()

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
  }).story(props =>
    `A list of Highly Connected Genes was obtained from the CTD output.`
  ).build()

export const Guilty_By_Association_Genes = MetaNode('Guilty_By_Association_Genes')
  .meta({
    label: `Extract Guilty By Association Genes`,
    description: 'Extract nodes that are "guilty by association" and connect your initial genes of interest within the graph.'
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
  }).story(props =>
    `A list of Guilty By Association Genes was obtained from the CTD output.`
  ).build()

export const CTD_Graph_Nodes = MetaNode('CTD_Graph_Nodes')
  .meta({
    label: 'CTD Graph',
    description: 'This is a visual display of the nodes that are found to be highly connected with CTD. Nodes that connect these nodes and are "guilty by association" will also be displayed.',
  })
  .inputs({ ctdResponseInfo: CTDResponseInfo })
  .output(GraphPlot)
  .resolve(async (props) => {
    let jsonGraph = props.inputs.ctdResponseInfo.jsonGraph;

    if(jsonGraph == null){
      throw new Error("No Gene Graph Nodes available, please use a different input gene set!");
    }

    return {
      edges: props.inputs.ctdResponseInfo.jsonGraph.interactions || [],
      nodes: props.inputs.ctdResponseInfo.jsonGraph.nodes.map(({ name, ...rest }) => ({ ...rest, id: name || '(no name)' })),
    }
  }).story(props =>
    `Graph Nodes were extracted from the CTD output.`
  ).build()

  export const GeneSet_CTD_String = MetaNode('GeneSet_CTD_String')
  .meta({
    label: `Connect the Dots in STRING`,
    description: 'Use CTD to "Connect the Dots" and identify highly connected set of proteins in the STRING protein interaction graph. *Please note 10-150 genes of interest are required to run CTD'
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
  }).story(props =>
    `CTD is applied which diffuses through all nodes in STRING [\\ref{doi:10.1093/nar/gku1003}] to identify nodes that are "guilty by association" and highly connected to the initial gene set of interest [\\ref{doi:10.1371/journal.pcbi.1009551}, \\ref{doi:10.1016/j.isci.2022.105799}].`
  ).build()

export const GeneSet_CTD_Wikipathways = MetaNode('GeneSet_CTD_Wikipathways')
  .meta({
    label: `Connect the Dots in Wikipathways`,
    description: 'Use CTD to "Connect the Dots" and identify highly connected set of genes in the WikiPathways pathway annotation graph. *Please note 10-150 genes of interest are required to run CTD'
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
  }).story(props =>
    `CTD is applied which diffuses through all nodes in WikiPathways [\\ref{doi:10.1093/nar/gkad960}] to identify nodes that are "guilty by association" and highly connected to the initial gene set of interest [\\ref{doi:10.1371/journal.pcbi.1009551}, \\ref{doi:10.1016/j.isci.2022.105799}].`
  ).build()

export const GenesFile_CTD_String = MetaNode('GenesFile_CTD_String')
  .meta({
    label: `CTD String For Genes Set File`,
    description: "Ensure a file contains a gene set, values separated by a \\n character  and with the extension .csv",
    icon: [file_transfer_icon]
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
  .story(props => ``)
  .build()

  export const GenesFile_CTD_Wikipathways = MetaNode('GenesFile_CTD_Wikipathways')
  .meta({
    label: `CTD Wikipathways For Genes Set File`,
    description: "Ensure a file contains a gene set, values separated by a \\n character  and with the extension .csv",
    icon: [file_transfer_icon]
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
  .story(props => ``)
  .build()
