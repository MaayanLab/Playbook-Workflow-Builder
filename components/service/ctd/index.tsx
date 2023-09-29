import { MetaNode } from '@/spec/metanode'
import { MultiFileURL } from '@/components/core/file'
import { GeneSet } from '@/components/core/input/set'
import { GraphPlot } from '@/components/viz/graph'
import { z } from 'zod'
import { FileURL } from '@/components/core/file'
import { fileAsStream } from  '@/components/core/file/api/download'
import { fileFromStream } from  '@/components/core/file/api/upload'
import { clientUploadFile } from  '@/components/core/file/api/upload/client'
import { downloadUrl } from '@/utils/download'
import { datafile_icon } from '@/icons'
import { Table, Cell, Column} from '@/app/components/Table'
import dynamic from 'next/dynamic'
import FormData from 'form-data'


const Matrix = dynamic(() => import('@/app/components/Matrix'))

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

//@ts-ignore
export async function getCTDFileResponse(formData): Promise<CTDResponse> {
  const res = await fetch(`http://genboree.org/pb-ctd/rest/playbook_ctd/ctd/file`, {
    method: 'POST',
    body: formData
  })
  return await res.json()
}

//@ts-ignore
export async function getCTDPrecalculationsResponse(formData): Promise<Readable> {
  return await fetch(`http://genboree.org/pb-ctd/rest/playbook_ctd/getCtdCustomMatrix`, {
    method: 'POST',
    body: formData
  });
}

export const CTD_Precalculations_RData = MetaNode('CTD_Precalculations_RData')
  .meta({
    label: 'CTD_Precalculations_RData',
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
  .inputs({ multiFileURL: MultiFileURL })
  .output(CTD_Precalculations_RData)
  .resolve(async (props) => {
    let geneListFile = props.inputs.multiFileURL[0];
    let adjMatrixFile = props.inputs.multiFileURL[1];

    const geneListFileReader: any = await fileAsStream(geneListFile);
    const adjMatrixFileReader: any = await fileAsStream(adjMatrixFile);
    const formData = new FormData();
    formData.append('geneList', geneListFileReader, geneListFile.filename);
    formData.append('customMatrix', adjMatrixFileReader, adjMatrixFile.filename);

    const respone = await getCTDPrecalculationsResponse(formData);
    if (!respone.ok) throw new Error(await respone.text)
    if (respone.body === null) throw new Error('An unexpected error occurred')
    //create a new file object from a stream
    const file = await fileFromStream(respone.body, `derived.${"customRDataFile.RData"}`)
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

export const Highly_Connected_Genes = MetaNode('Highly_Connected_Genes')
  .meta({
    label: `Highly_Connected_Genes`,
    description: "Get a list of Highly_Connected_Genes from the CTD response."
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
    `Get a list of Highly_Connected_Genes from the CTD response.`
  ).build()


export const Guilty_By_Association_Genes = MetaNode('Guilty_By_Association_Genes')
  .meta({
    label: `Guilty_By_Association_Genes`,
    description: "Get a list of Guilty_By_Association_Genes from the CTD response."
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
    `Get a list of Guilty_By_Association_Genes from the CTD response.`
  ).build()


export const CTD_Graph_Nodes = MetaNode('CTD_Graph_Nodes')
  .meta({
    label: `CTD_Graph_Nodes`,
    description: "CTD_Graph_Nodes."
  })
  .inputs({ ctdResponseInfo: CTDResponseInfo })
  .output(GraphPlot)
  .resolve(async (props) => {
    let jsonGraphObj = props.inputs.ctdResponseInfo.jsonGraph;

    jsonGraphObj.edges = jsonGraphObj.interactions;
    delete jsonGraphObj.interactions;

    for(let n in jsonGraphObj.nodes){
      let node = jsonGraphObj.nodes[n];
      node.id = node.name;
      delete node.name;
    }
    
    return jsonGraphObj;
  }).story(props =>
    `CTD_Graph_Nodes.`
  ).build()
