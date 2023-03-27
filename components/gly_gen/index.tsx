import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { GeneInfo } from '../service/mygeneinfo'
import { z } from 'zod'
import { glygen_icon } from '@/icons'

export const GlyGenResponse = z.object({
  queryinfo: z.object({
    query: z.object({
      recommended_gene_name: z.string()
    }),
  }),
  results: z.array(z.object({
    uniprot_canonical_ac: z.string()
  }))
})

export type GlyGenResponseType = z.infer<typeof GlyGenResponse>

export const GlyGenResponseNode = MetaNode('GlyGenResponse')
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

export const ProteinProductInformation = MetaNode('ProteinProductInformation')
  .meta({
    label: 'Search Glygen for Protein Products',
    description: 'Protein product records in GlyGen for the gene',
    icon: [glygen_icon],
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
    const response = GlyGenResponse.parse(await request.json())
    return response
  })
  .story(props =>
    `Next, the GlyGen database [\\ref{York, William S et al. “GlyGen: Computational and Informatics Resources for Glycoscience.” Glycobiology vol. 30,2 (2020): 72-73. doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from the gene.`
  )
  .build()
