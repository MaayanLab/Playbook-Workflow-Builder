import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { GeneSymbol } from '@/components/Gene'

export type MyGeneInfoHit = {
  _id: string,
  _score: number,
  entrezgene: string,
  name: string,
  symbol?: string,
  taxid: number,
}

export type MyGeneInfo = {
  _id: string,
  entrezgene: string,
  name: string,
  symbol: string,
  taxid: number,
}

async function mygeneinfo_query(geneId: string): Promise<{total: number, hits: Array<MyGeneInfoHit>}> {
  const res = await fetch(`https://mygene.info/v3/query?q=${encodeURIComponent(geneId)}`)
  return await res.json()
}

async function mygeneinfo(geneId: string): Promise<MyGeneInfo> {
  const res = await fetch(`https://mygene.info/v3/gene/${encodeURIComponent(geneId)}`)
  return await res.json()
}

export const GeneInfo = MetaNode.createData('GeneInfo')
  .meta({
    label: 'Gene Information',
    description: 'A Gene resolved with MyGeneInfo',
  })
  .codec<MyGeneInfo>()
  .view(geneinfo => (
    <div>
      <a href={`https://www.ncbi.nlm.nih.gov/gene/${geneinfo.entrezgene}`}>{geneinfo.symbol}</a> {geneinfo.name}
    </div>
  ))
  .build()

export const GeneInfoFromGeneSymbol = MetaNode.createProcess('GeneInfoFromGeneSymbol')
  .meta({
    label: 'Resolve Gene from Symbol',
    description: 'Resolve gene info from gene symbol with MyGeneInfo',
  })
  .codec()
  .inputs({ gene: GeneSymbol })
  .output(GeneInfo)
  .resolve(async (props) => {
    const results = await mygeneinfo_query(props.inputs.gene)
    const hits = results.hits.filter(hit => hit.symbol)
    const exactMatch = hits.filter(hit => hit.symbol.toUpperCase() == props.inputs.gene.toUpperCase())[0]
    const _id: string | undefined = exactMatch !== undefined ? exactMatch._id : hits[0]._id
    if (_id === undefined) throw new Error(`Could not identify a gene for the symbol ${props.inputs.gene} in mygene.info`)
    return await mygeneinfo(_id)
  })
  .build()
