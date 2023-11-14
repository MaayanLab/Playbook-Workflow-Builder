import { MetaNode } from '@/spec/metanode'
import { FileURL } from '@/components/core/file'
import { CTDPrecalculationsFileURLs, CTDUseCustomMatrixFileURLs } from './input'
import { GeneSet } from '@/components/core/input/set'
import { z } from 'zod'
import { file_transfer_icon, datafile_icon } from '@/icons'
import { fileAsStream } from  '@/components/core/file/api/download'
import { GraphPlot } from '@/components/viz/graph'
import { fileFromStream } from  '@/components/core/file/api/upload'
import { downloadUrl } from '@/utils/download'
import { Table, Cell, Column} from '@/app/components/Table'
import FormData from 'form-data'
import axios from 'axios'
import { Readable } from 'stream'

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
  const res = await axios.post(`http://genboree.org/pb-ctd/rest/playbook_ctd/ctd/file`, formData, {
    headers: { ...formData.getHeaders() },
    responseType: 'json',
  })
  return res.data
}

export async function getCTDPrecalculationsResponse(formData: FormData): Promise<Readable> {
  const res = await axios.post(`http://genboree.org/pb-ctd/rest/playbook_ctd/getCtdCustomMatrix`, formData, {
    headers: { ...formData.getHeaders() },
    responseType: 'stream',
  });
  if (res.status !== 200) throw new Error(await res.data.read())
  return res.data
}

export async function getCTDUseCustomMatrix(formData: FormData): Promise<CTDResponse> {
  const res = await axios.post(`http://genboree.org/pb-ctd/rest/playbook_ctd//ctd/useCustomMatrix`, formData, {
    headers: { ...formData.getHeaders() },
    responseType: 'json',
  });
  return res.data
}

export const CTD_FileDownload = MetaNode('CTD_FileDownload')
  .meta({
    label: 'CTD_FileDownload',
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
    label: `Execute_CTD_Precalculations`,
    description: "Execute CTD Precalculations using a custom geme list and ajd. matrix file in order to get am RData file for the final step."
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

    const respone = await getCTDPrecalculationsResponse(formData);
    //create a new file object from a stream
    const file = await fileFromStream(respone, `derived.${"customRDataFile.RData"}`)
    return file;
  }).story(props =>   
    "The two files where send to the CTD API for precalculations."
  ).build()

  export const CTDResponseInfo = MetaNode('CTDResponseInfo')
  .meta({
    label: 'CTDResponseInfo',
    description: '',
    icon: [datafile_icon]
  })
  .codec(CTDResponseC)
  .view(props => {
    return(
      <div>
        <p>Highly Connected Genes num: {props.highlyConnectedGenes.length}</p>
        <p>Guilty By Association Genes num: {props.guiltyByAssociationGenes.length}</p>
        <p>Json Graph Nodes num: {props.jsonGraph.nodes.length}</p>
      </div>
    )
  }).build()

  export const CTD_UseCustomMatrix = MetaNode('CTD_UseCustomMatrix')
  .meta({
    label: `CTD_UseCustomMatrix`,
    description: "Execute CTD Precalculations using a custom geme list and ajd. matrix file in order to get am RData file for the final step."
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

    return await getCTDUseCustomMatrix(formData);
  }).story(props =>   
    "The three files where send to the CTD API for precalculations."
  ).build()

export const Highly_Connected_Genes = MetaNode('Highly_Connected_Genes')
  .meta({
    label: `Extract Highly Connected Genes`,
    description: "Get a list of Highly Connected Genes from the CTD output"
  })
  .inputs({ ctdResponseInfo: CTDResponseInfo })
  .output(GeneSet)
  .resolve(async (props) => {
    let geneSet = {
      description: 'Highly_Connected_Genes, CTD response',
      set:  props.inputs.ctdResponseInfo.highlyConnectedGenes
    };
    return geneSet;
  }).story(props =>
    `A list of Highly Connected Genes was obtained from the CTD output.`
  ).build()

export const Guilty_By_Association_Genes = MetaNode('Guilty_By_Association_Genes')
  .meta({
    label: `Extract Guilty By Association Genes`,
    description: "Get a list of Guilty By Association Genes from the CTD output"
  })
  .inputs({ ctdResponseInfo: CTDResponseInfo })
  .output(GeneSet)
  .resolve(async (props) => {
    let geneSet = {
      description: 'Guilty_By_Association_Genes, CTD response',
      set:  props.inputs.ctdResponseInfo.guiltyByAssociationGenes
    };
    return geneSet;
  }).story(props =>
    `A list of Guilty By Association Genes was obtained from the CTD output.`
  ).build()

export const CTD_Graph_Nodes = MetaNode('CTD_Graph_Nodes')
  .meta({
    label: 'CTD Graph',
    description: 'A graph showing the CTD output.',
  })
  .inputs({ ctdResponseInfo: CTDResponseInfo })
  .output(GraphPlot)
  .resolve(async (props) => {
    return {
      edges: props.inputs.ctdResponseInfo.jsonGraph.interactions || [],
      nodes: props.inputs.ctdResponseInfo.jsonGraph.nodes.map(({ name, ...rest }) => ({ ...rest, id: name || '(no name)' })),
    }
  }).story(props =>
    `Graph Nodes were extracted from the CTD output.`
  ).build()

/* kegg is not used for now because of licencing issues
export const GeneSet_CTD_Kegg = MetaNode('GeneSet_CTD_Kegg')
  .meta({
    label: `GeneSet_CTD_Kegg`,
    description: "Get a CTD response for a set of genes for graph type kegg."
  })
  .inputs({ geneset: GeneSet })
  .output(CTDResponseInfo)
  .resolve(async (props) => {
    let requestBody = {
      "graphType": "kegg",
      "geneList": props.inputs.geneset.set
    }
    return await getCTDGenSetResponse(JSON.stringify(requestBody));
  }).story(props =>
    `Get a CTD response for a set of genes for graph type kegg.`
  ).build()
*/

export const GeneSet_CTD_String = MetaNode('GeneSet_CTD_String')
  .meta({
    label: `GeneSet_CTD_String`,
    description: "Get a CTD response for a set of genes for graph type string."
  })
  .inputs({ geneset: GeneSet })
  .output(CTDResponseInfo)
  .resolve(async (props) => {
    let requestBody = {
      "graphType": "string",
      "geneList": props.inputs.geneset.set
    }
    return await getCTDGenSetResponse(JSON.stringify(requestBody));
  }).story(props =>
    `Get a CTD response for a set of genes for graph type string.`
  ).build()

/* //kegg is not used for now because of licencing issues
export const GenesFile_CTD_Kegg = MetaNode('GenesFile_CTD_Kegg')
  .meta({
    label: `GenesFile_CTD_Kegg`,
    description: "Ensure a file contains a gene set, values separated by a \\n character  and with the extension .csv",
    icon: [file_transfer_icon]
  })
  .inputs({ file: FileURL })
  .output(CTDResponseInfo)
  .resolve(async (props) => {
    const fileReader = await fileAsStream(props.inputs.file);

    const formData = new FormData();
    formData.append('csvGenesFile', fileReader, props.inputs.file.filename);
    formData.append('graphType', "kegg");
    return await getCTDFileResponse(formData);
  }).build()*/

export const GenesFile_CTD_String = MetaNode('GenesFile_CTD_String')
  .meta({
    label: `GenesFile_CTD_String`,
    description: "Ensure a file contains a gene set, values separated by a \\n character  and with the extension .csv",
    icon: [file_transfer_icon]
  })
  .inputs({ file: FileURL })
  .output(CTDResponseInfo)
  .resolve(async (props) => {
    const fileReader = await fileAsStream(props.inputs.file);

    const formData = new FormData();
    formData.append('csvGenesFile', fileReader, props.inputs.file.filename);
    formData.append('graphType', "string");
    return await getCTDFileResponse(formData);
  })
  .story(props => ``)
  .build()
