import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { GeneInfo } from '../service/mygeneinfo'
import { z } from 'zod'

export const GlyGenGeneType = {
  gene_name: typeof GeneInfo
}

export const GlyGenResponse = z.object({
  queryinfo: z.any(),
  results: z.array(z.any())
})

export type GlyGenResponseType = z.infer<typeof GlyGenResponse>

export const GlyGenResponseNode = MetaNode.createData('GlyGenResponse')
  .meta({
    label: 'GlyGen Response',
    description: 'GlyGen response object',
  })
  .codec(GlyGenResponse)
  .view(data => (
    <div>
      Query:
      <pre >
        {JSON.stringify(data.queryinfo.query, undefined, 2)}
      </pre>
      Results:
      {data.results.map((result, index) => (
        <div key={index}>
          {JSON.stringify(result.uniprot_canonical_ac, undefined, 2)}
        </div>
      ))}
    </div>
  ))
  .build()

export const ProteinProductInformation = MetaNode.createProcess('ProteinProductInformation')
  .meta({
    label: 'Protein Product Information',
    description: 'Search for protein records in GlyGen',
  })
  .inputs({ gene: GeneInfo })
  .output(GlyGenResponseNode)
  .resolve(async (props) => {
    const query = encodeURIComponent(JSON.stringify({
      recommended_gene_name: props.inputs.gene.symbol,
    }))
    const request = await fetch(`https://api.glygen.org/directsearch/protein/?query=${query}`, {
      method: 'GET',
      headers: {
        accept: 'application/json'
      }
      })
    const response = await request.json()
    return response
  })
  .build()