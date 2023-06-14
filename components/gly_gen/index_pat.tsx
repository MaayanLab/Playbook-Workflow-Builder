import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { GeneInfo, GeneInfoFromGeneTerm } from '../service/mygeneinfo'
import { z } from 'zod'
import { glygen_icon } from '@/icons'
import { GeneTerm } from '@/components/core/input/term'

export const GlyGenIDtoResponse = z.object({
  list_id: z.string(),
})

export const GlyGenResponse = z.object({ })
//   queryinfo: z.object({
//     query: z.object({
//       gene_name: z.string()
//     }),
//   }),
//   results: z.array(z.object({
//     uniprot_canonical_ac: z.string()
//   }))
// })

export type GlyGenResponseType = z.infer<typeof GlyGenResponse>
export type GlyGenIDType = z.infer<typeof GlyGenIDtoResponse>

export const GlyGenResponseNode = MetaNode('GlyGenResponse')
  .meta({
    label: 'GlyGen Protein Products',
    description: 'Protein product records in GlyGen',
    icon: [glygen_icon],
  })
  .codec(GlyGenResponse)
  .view(data => (
    <div>
      Query:
      <pre >
        {/* {JSON.stringify(data.queryinfo.query, undefined, 2)} */}
        {JSON.stringify(data, undefined, 2)}
      </pre>
      {/* Results:
      {data.results.map((result, index) => (
        <div key={index}>
          {JSON.stringify(result.uniprot_canonical_ac, undefined, 2)}
        </div> */}
      {/* ))} */}
    </div>
  ))
  .build()

async function SearchGlyGenProteins(gene_name: string) {
  // TODO: Breakout ID request here
  const id_request = await fetch('https://api.glygen.org/directsearch/protein/search', {
    method: 'POST',
    headers: {
      accept: 'application/json',
      'Content-Type': 'application/json',
    },
    body: JSON.stringify({ gene_name: gene_name }),
    })
  return id_request.json()
}

async function ResolveGlyGenProteinByID(gene_id: string){
  // TODO: Breakout retrieval of info here
  const protein_response = await fetch('https://api.glygen.org/directsearch/protein/list', {
    method: 'POST',
    headers: {
      accept: 'application/json',
      'Content-Type': 'application/json',
    },
    body: JSON.stringify({ id: gene_id }),
    })
  return protein_response
}

async function SearchAndResolve(gene_name: string) {
  try {
    const glygen_id = await SearchGlyGenProteins(gene_name);
    const response_body = await ResolveGlyGenProteinByID(glygen_id['list_id']);
    return response_body.json()
  } catch(e) {
    console.log(e)
    throw e
  }
}

export const ProteinProductInformation = MetaNode('ProteinProductInformation')
  .meta({
    label: 'Search Glygen for Protein Products',
    description: 'Find protein product records in GlyGen for the gene',
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ gene: GeneInfo })
  .output(GlyGenResponseNode)
  .resolve(async (props) => {
    // const query = encodeURIComponent(JSON.stringify({
    //   gene_name: props.inputs.gene.symbol,
    // }))
    // const id_request = await fetch('https://api.glygen.org/directsearch/protein/search', {
    //   method: 'POST',
    //   headers: {
    //     accept: 'application/json',
    //     'Content-Type': 'application/json',
    //   },
    //   body: JSON.stringify({ gene_name: props.inputs.gene.symbol }),
    //   })
    // const protein_response = await fetch('https://api.glygen.org/directsearch/protein/list', {
    //   method: 'POST',
    //   headers: {
    //     accept: 'application/json',
    //     'Content-Type': 'application/json',
    //   },
    //   body: JSON.stringify({ id: GlyGenIDtoResponse.parse(await id_request.json())['list_id'] }),
    //   })
    await SearchAndResolve(props.inputs.gene.symbol).then(response => {
      return response
    })
    // const response = GlyGenResponse.parse()
    // return response
  })
  .story(props =>
    `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from ${props.inputs.gene.symbol}.`
  )
  .build()

export const ProteinProductInformationFromGene = MetaNode('ProteinProductInformationFromGene')
  .meta(ProteinProductInformation.meta)
  .inputs({ gene: GeneTerm })
  .output(ProteinProductInformation.output)
  .resolve(async (props) => {
    const gene = await GeneInfoFromGeneTerm.resolve(props)
    return await ProteinProductInformation.resolve({ inputs: { gene } })
  })
  .story(props =>
    `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from ${props.inputs.gene}.`
  )
  .build()
