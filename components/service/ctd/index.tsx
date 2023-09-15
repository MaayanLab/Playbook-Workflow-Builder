import { MetaNode } from '@/spec/metanode'
import { FileURL } from '@/components/core/file'
import { GeneSet } from '@/components/core/input/set'
import { z } from 'zod'
import { file_transfer_icon, datafile_icon } from '@/icons'
import { fileAsStream } from  '@/components/core/file/api/download'
import { GraphPlot } from '@/components/viz/graph'
import FormData from 'form-data'

const CTDResponseInfoC = z.object({
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
})
export type CTDResponseInfo = z.infer<typeof CTDResponseInfoC>

export async function getCTDFileResponse(formData: any): Promise<CTDResponseInfo> {
  const res = await fetch(`http://genboree.org/pb-ctd/rest/playbook_ctd/ctd/file`, {
    method: 'POST',
    body: formData
  })
  return await res.json()
}

export async function getCTDGenSetResponse(strValue: string): Promise<CTDResponseInfo> {
  const res = await fetch(`http://genboree.org/pb-ctd/rest/playbook_ctd/ctd/geneList`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
    },
    body: strValue
  })
  return await res.json()
}

export const CTDResponseInfo = MetaNode('CTDResponseInfo')
  .meta({
    label: 'CTD Results',
    description: 'The results of CTD',
    icon: [datafile_icon]
  })
  .codec(CTDResponseInfoC)
  .view(props => {
    return(
      <div>
        <p>Highly Connected Genes num: {props.highlyConnectedGenes.length}</p>
        <p>Guilty By Association Genes num: {props.guiltyByAssociationGenes.length}</p>
        <p>Json Graph Nodes num: {props.jsonGraph.nodes.length}</p>
      </div>
    )
  }).build()

/* kegg is not used for now because of licencing issues
export const GeneSet_CTD_Kegg = MetaNode('GeneSet_CTD_Kegg')
  .meta({
    label: `Perform CTD with Gene Set using KEGG`,
    description: "Get a CTD response for a set of genes using KEGG"
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
    `CTD was performed with the gene set using KEGG.`
  ).build() */


export const GeneSet_CTD_String = MetaNode('GeneSet_CTD_String')
  .meta({
    label: `Perform CTD with Gene Set using STRING`,
    description: "Get a CTD response for a set of genes using STRING"
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
    `CTD was performed with the gene set using STRING.`
  ).build()


// TODO: maybe a Gene Set From File can replace this?
/* //kegg is not used for now because of licencing issues
export const GenesFile_CTD_Kegg = MetaNode('GenesFile_CTD_Kegg')
  .meta({
    label: `Perform CTD with File using KEGG`,
    description: 'File should contain a single column gene set in csv format',
    icon: [file_transfer_icon]
  })
  .inputs({ file: FileURL })
  .output(CTDResponseInfo)
  .resolve(async (props) => {
    const fileReader: any = await fileAsStream(props.inputs.file);

    const formData = new FormData();
    formData.append('csvGenesFile', fileReader, props.inputs.file.filename);
    formData.append('graphType', "kegg");
    return await getCTDFileResponse(formData);
  }).story(props =>
    `CTD was performed with the gene set using KEGG.`
  ).build()*/

export const GenesFile_CTD_String = MetaNode('GenesFile_CTD_String')
  .meta({
    label: `Perform CTD with File using STRING`,
    description: 'File should contain a single column gene set in csv format.',
    icon: [file_transfer_icon]
  })
  .inputs({ file: FileURL })
  .output(CTDResponseInfo)
  .resolve(async (props) => {
    const fileReader: any = await fileAsStream(props.inputs.file);

    const formData = new FormData();
    formData.append('csvGenesFile', fileReader, props.inputs.file.filename);
    formData.append('graphType', "string");
    return await getCTDFileResponse(formData);
  }).story(props =>
    `CTD was performed with the gene set using STRING.`
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

//just for demo, this needs to be presented with a proper graph visualization
export const CTDGraph = MetaNode('CTDGraph')
  .meta({
    label: 'CTD Graph',
    description: 'A graph showing the CTD output.',
    icon: [datafile_icon]
  })
  .codec()
  .view(props => {
    return(
      <div>
          {JSON.stringify(props)}
      </div>
    )
  }).build()

export const CTD_Graph_Nodes = MetaNode('CTD_Graph_Nodes')
  .meta({
    label: `Extract the CTD Graph Nodes`,
    description: "CTD Graph Nodes"
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
