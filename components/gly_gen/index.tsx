import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { GeneInfo, GeneInfoFromGeneTerm } from '../service/mygeneinfo'
import { z } from 'zod'
import { glygen_icon, protein_icon } from '@/icons'
import { GeneTerm } from '@/components/core/input/term'
import { filterGlyGenResults } from './utils'
import { Properties } from '@blueprintjs/icons/lib/esm/generated/16px/paths'
import { filter } from '@/utils/dict'

export const GlyGenResponse = z.object({
    gene: z.object({
      name: z.string(),
      locus: z.object({
        chromosome: z.string(),
        start_pos: z.number(),
        end_pos: z.number(),
        strand: z.string()
      })
    }),
    uniprot: z.object({
      uniprot_id: z.string(),
      uniprot_canonical_ac: z.string(),
      length: z.number()
    }),
    protein_names: z.object({
      name: z.string()
    }),
    // mass: z.object({
    //   chemical_mass: z.string()
    // }),
    // refseq_ac: z.string(),
    // refseq_name: z.string(),
    species: z.object({
      name: z.string(),
      common_name: z.string(),
      taxid: z.string(),
    })
  })

// export const GlycanResponse = z.object({})
// TODO - List API
// export const GlyGenListResponse = z.object({
//   results: z.array(z.object({
//     gene_name: z.string(),
//     protein_name: z.string(),
//     uniprot_canonical_ac: z.string()
//   }))
// })
export const ProteinResponse = z.object({})

export type GlyGenResponseType = z.infer<typeof GlyGenResponse>
export type ProtienResponseType = z.infer<typeof ProteinResponse>

export const ProteinResponseNode = MetaNode('Protein')
  .meta({
    label: 'Basic Protein',
    description: 'Protein primitive type sourced from the Glygen server',
    icon: [protein_icon],
  })
  .codec(ProteinResponse)
  .view( () => (
    <div>
      Results: Pending
    </div>)
  )
  .build()

export const GlyGenResponseNode = MetaNode('GlyGenResponse')
  .meta({
    label: 'GlyGen Protein Products',
    description: 'Protein product records in GlyGen',
    icon: [glygen_icon],
  })
  .codec(GlyGenResponse)
  .view( data => (
    <div>
        <div>Gene Name: <b>{data.gene.name}</b></div>
        <div>Gene location: Chromosome: {data.gene.locus.chromosome} ({data.gene.locus.start_pos} - {data.gene.locus.end_pos}, '{data.gene.locus.strand}' strand)</div>
        <div>UniProtKB ID: {data.uniprot.uniprot_id}</div>
        <div>UniProtKB Accession: {data.uniprot.uniprot_canonical_ac}</div>
        <div>Protein Length: <b>{data.uniprot.length}</b></div>
        <div>UniProtKB Protein Name(s): {data.protein_names.name}</div>
        <div>Organism: <b>{data.species.name} ({data.species.common_name}; TaxID: {data.species.taxid})</b></div>
    </div>
  ))
  // .view( data => {
    // <div>
    //   Results:
    //   {data.results.map((result, index) => (
    //     <div>
    //         <div>Gene Name: {result.gene_name}</div>
    //         <div>Protien Name: {result.protein_name}</div>
    //         <div>UniProtKB Accession: {result.uniprot_canonical_ac}</div>
    //     <br/>
    //     </div>
    //   ))}
    // </div>
    // }
  // ))
  .build()

// function filterGlyGenResults(result, gene_name) {
//   for (const idx in result['results']) {
//     if (result['results'][idx].hasOwnProperty('gene_name')){
//       const geneName = result['results'][idx]['gene_name'].toUpperCase();
//       const speciesName = result['results'][idx]['organism'];
//       if (geneName === gene_name.toUpperCase() && speciesName === 'Homo sapiens'){ 
//           console.log('==========================');
//           console.log(`Human result: ${result['results'][idx]['gene_name']}`)
//           console.log('==========================');
//           return result['results'][idx]
//       }
//     }
//   }
// }

export const ProteinProductInformation = MetaNode('ProteinProductInformation')
  .meta({
    label: 'Search Glygen for Protein Products',
    description: 'Find protein product records in GlyGen for the gene',
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ gene: GeneInfo })
  .output(ProteinResponseNode)
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
    const searchResult = await protein_response.json()
    const filteredResult = filterGlyGenResults(searchResult, props.inputs.gene.symbol);
    return filteredResult;
  })
  // .prompt(() => {

  // })
  .story(props =>
    `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from ${props.inputs ? props.inputs.gene.symbol : 'the gene'}.`
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
    `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from ${props.inputs ? props.inputs.gene : 'the gene'}.`
  )
  .build()
