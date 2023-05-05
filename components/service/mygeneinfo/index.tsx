import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { FileURL } from '@/components/core/file'
import { GeneTerm } from '@/components/core/input/term'
import { RegulatoryElementSet } from '@/components/core/input/set'
import { z } from 'zod'
import { gene_icon, mygeneinfo_icon, file_transfer_icon, datafile_icon } from '@/icons'

export const MyGeneInfoHitC = z.object({
  hits: z.array(
    z.object({
      _id: z.string(),
      _score: z.number(),
      symbol: z.string().optional(),
      name: z.string().optional(),
      taxid: z.number().optional(),
      entrezgene: z.string().optional(),
    }),
  ),
})
export type MyGeneInfoHit = z.infer<typeof MyGeneInfoHitC>

export const MyGeneInfoC = z.object({
  _id: z.string(),
  symbol: z.string(),
  entrezgene: z.string().optional(),
  name: z.string().optional(),
  taxid: z.number().optional(),
  ensembl: z.object({
    gene: z.string(),
  }).optional()
})
export type MyGeneInfo = z.infer<typeof MyGeneInfoC>

export const MyGeneInfoByTermC = z.object({
  data: z.object({
    ld: z.object({
      RegulatoryElement: z.array(z.object({ 
            entId: z.string(),
            ldhId: z.string()
          })
        )
    })
  })
})
export type MyGeneInfoByTerm = z.infer<typeof MyGeneInfoByTermC>

async function mygeneinfo_query(geneId: string): Promise<{total: number, hits: Array<MyGeneInfoHit>}> {
  const res = await fetch(`https://mygene.info/v3/query?q=${encodeURIComponent(geneId)}`)
  return await res.json()
}

async function mygeneinfo(geneId: string): Promise<MyGeneInfo> {
  const res = await fetch(`https://mygene.info/v3/gene/${encodeURIComponent(geneId)}`)
  return await res.json()
}

async function myGeneInfoByGeneTerm(geneTerm: string): Promise<MyGeneInfoByTerm> {
  const res = await fetch(`https://genboree.org/cfde-gene-dev/Gene/id/${encodeURIComponent(geneTerm)}`)
  return await res.json()
}

export const GeneInfo = MetaNode('GeneInfo')
  .meta({
    label: 'Gene Information',
    description: 'A Gene resolved with MyGeneInfo',
    icon: [gene_icon],
  })
  .codec(MyGeneInfoC)
  .view(geneinfo => (
    <div>
      <a href={`https://www.ncbi.nlm.nih.gov/gene/${geneinfo.entrezgene}`}>{geneinfo.symbol}</a> {geneinfo.name}
    </div>
  ))
  .build()

export const GeneInfoFromGeneTerm = MetaNode('GeneInfoFromGeneTerm')
  .meta({
    label: 'Resolve Gene Info from Term',
    description: 'Resolve gene info from gene term with MyGeneInfo',
    icon: [mygeneinfo_icon],
  })
  .inputs({ gene: GeneTerm })
  .output(GeneInfo)
  .resolve(async (props) => {
    return await getGeneData(props.inputs.gene);
  })
  .build()

  export async function getGeneData(geneSymbol: string){
    const results = MyGeneInfoHitC.parse(await mygeneinfo_query(geneSymbol))
    const hits = results.hits.filter((hit): hit is MyGeneInfoHit['hits'][0] & { symbol: string } => !!hit.symbol)
    const exactMatch = hits.filter(hit => hit.symbol.toUpperCase() == geneSymbol.toUpperCase())[0]
    const _id: string | undefined = exactMatch !== undefined ? exactMatch._id : hits[0]._id
    if (_id === undefined) throw new Error(`Could not identify a gene for the symbol ${geneSymbol} in mygene.info`)
    return await mygeneinfo(_id)
  }


  export const GetRegulatoryElementsForGeneInfo = MetaNode('GetRegulatoryElementsForGeneInfo')
  .meta({
    label: 'Resolve Reg. Elements from Gene Info',
    description: 'GetRegulatoryElementsForGeneInfo',
  })
  .inputs({ geneInfo: GeneInfo  })
  .output(RegulatoryElementSet)
  .resolve(async (props) => {
    const response =  await myGeneInfoByGeneTerm(props.inputs.geneInfo.symbol);
    if(response.data == null || response.data.ld == null){
      return [];
    }
    return response.data.ld.RegulatoryElement.map(({ entId }) => entId );
  })
  .build()
  
  const CTDResponseInfoC = z.object({ 
      "highlyConnectedGenes": z.any(),
      "guiltyByAssociationGenes": z.any(),
      "jsonGraph": z.object({ 
        "nodes": z.array(z.object({
          "name": z.string(),
          "type": z.string()
        })),
        "interactions": z.array(z.object({
          "source": z.string(),
          "target": z.string()
        })),
      })
  });
  export type CTDResponseInfo = z.infer<typeof CTDResponseInfoC>

  
  export async function getCTDResponse(formData: FormData): Promise<CTDResponseInfo> {
    const res = await fetch(`http://5.161.50.225:8018/rest/playbook_ctd/ctd/file`, {
      method: 'POST',
      headers: {
        'Content-Type': 'multipart/form-data',
      },
      body: formData
    })
    return await res.json()
  }

  export const CTDResponseInfo = MetaNode('CTDResponseInfo')
  .meta({
    label: 'CTDResponseInfo',
    description: '',
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
/*
        <p>Highly Connected Genes num: {props.highlyConnectedGenes.length}</p>
        <p>Guilty By Association Genes num: {props.guiltyByAssociationGenes.length}</p>
        <p>Json Graph Nodes num: {props.jsonGraph.nodes.length}</p>
*/






export const GeneCTD = MetaNode('GeneCTD')
  .meta({
    label: `GeneCTD`,
    description: "Ensure a file contains a gene set, values separated by a \\n character  and with the extension .csv",
    icon: [file_transfer_icon]
  })
  .inputs({ file: FileURL })
  .output(CTDResponseInfo)
  .resolve(async (props) => {
    const formData = new FormData();
 
    async function download(file: { url: string }) {
      if (file.url.startsWith('file://')) {
        const fs = typeof window === 'undefined' ? require('fs') : undefined;
        return fs.readFileSync(file.url.slice('file://'.length)).toString();
      } else {
        const req = await fetch(file.url);
        return await req.text();
      }
    }

    const myFileContent = await download({"url": props.inputs.file});
    //let file = new File([myFileContent], "CTD_test.csv", { type: "text/html" });
    formData.append('csvGenesFile', new Blob([myFileContent], { type: "text/html" }), "CTD_test.csv");

    return await getCTDResponse(formData);
  })
  .build()


  
  
   
