import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { GeneInfo, GeneInfoFromGeneTerm } from '../service/mygeneinfo'
import { z } from 'zod'
import { glygen_icon } from '@/icons'
import { GeneTerm } from '@/components/core/input/term'

export const GlyGenResponse = z.object({
    results: z.array(z.object({
        gene_name: z.string(),
        protein_name: z.string(),
        uniprot_canonical_ac: z.string()
    }))
})

export type GlyGenResponseType = z.infer<typeof GlyGenResponse>

export const GlyGenResponseNode = MetaNode('GlyGenResponse')
  .meta({
    label: 'GlyGen Protein Products',
    description: 'Protein product records in GlyGen',
  })
  .codec(GlyGenResponse)
  .view(data => (
    <div>
      Results:
      {data.results.map((result, index) => (
        <div>
            <div>Gene Name: {result.gene_name}</div>
            <div>Protien Name: {result.protein_name}</div>
            <div>UniProtKB Accession: {result.uniprot_canonical_ac}</div>
        <br/>
        </div>
      ))}
    </div>
  ))
  .build()

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
    const id_request = await fetch('https://api.glygen.org/protein/search/', {
      method: 'POST',
      headers: {
        accept: 'application/json',
       'Content-Type': 'application/json',
      },
      body: JSON.stringify({ gene_name: props.inputs.gene.symbol }),
    })
    const id = await id_request.json()
    const protein_response = await fetch('https://api.glygen.org/protein/list/', {
      method: 'POST',
      headers: {
        accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ id: id['list_id'] }),
      })
    return protein_response.json()
  })
  .story(props =>
    `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from ${props.inputs.gene.symbol}.`
  )
  .build()
export const ProteinProductInformationFromGene = MetaNode('ProteinProductInformationFromGene')
  .meta({
    label: 'Search Glygen for Protein Products',
    description: 'Find protein product records in GlyGen for the gene',
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ gene: GeneTerm })
  .output(GlyGenResponseNode)
  .resolve(async (props) => {
    const gene = await GeneInfoFromGeneTerm.resolve(props)
    return await ProteinProductInformation.resolve({ inputs: { gene } })
  })
  .story(props =>
    `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from ${props.inputs.gene}.`
  )
  .build()