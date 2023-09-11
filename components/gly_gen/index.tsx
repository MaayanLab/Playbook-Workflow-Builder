import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { GeneInfo, GeneInfoFromGeneTerm } from '../service/mygeneinfo'
import { z } from 'zod'
import { glygen_icon, protein_icon } from '@/icons'
import { GeneTerm, ProteinTerm } from '@/components/core/input/term'
import { filterGlyGenResults } from './utils'
import { Properties } from '@blueprintjs/icons/lib/esm/generated/16px/paths'
import { filter } from '@/utils/dict'

export const GlyGenProteinResponse = z.object({
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
    }),
    glycoprotein: z.object({
      glycosylation: z.boolean()
    })
  })

// TODO - List API(s)

export type GlyGenProteinResponseType = z.infer<typeof GlyGenProteinResponse>
// export type ProteinResponseType = z.infer<typeof ProteinResponse>

// export const ProteinResponseNode = MetaNode('Protein')
//   .meta({
//     label: 'Basic Protein',
//     description: 'Protein primitive type sourced from the GlyGen server',
//     icon: [protein_icon],
//   })
//   .codec(ProteinResponse)
//   .view( () => (
//     <div>
//       Results: Pending
//     </div>)
//   )
//   .build()

export const GlyGenProteinResponseNode = MetaNode('GlyGenProteinResponse')
  .meta({
    label: 'GlyGen Protein Products',
    description: 'Protein product records in GlyGen',
    icon: [glygen_icon],
  })
  .codec(GlyGenProteinResponse)
  .view( data => (
    <div>
        <div>Gene Name: <b>{data.gene.name}</b></div>
        <div>UniProtKB Accession: {data.uniprot.uniprot_canonical_ac}</div>
        <div>Gene location: Chromosome: {data.gene.locus.chromosome} ({data.gene.locus.start_pos} - {data.gene.locus.end_pos}, '{data.gene.locus.strand}' strand)</div>
        <div>UniProtKB ID: {data.uniprot.uniprot_id}</div>
        <div>Protein Length: <b>{data.uniprot.length}</b></div>
        <div>UniProtKB Protein Name(s): {data.protein_names.name}</div>
        <div>Organism: <b>{data.species.name} ({data.species.common_name}; TaxID: {data.species.taxid})</b></div>
        <div>Glycoprotein: {data.glycoprotein.glycosylation ? 'True' : 'False'}</div>
    </div>
  ))
  .build()

export const GlyGenProtein = MetaNode('GGP')
.meta({
  label: 'Search GlyGen by Protein Name for Protein Products',
  description: 'Find protein product records in GlyGen for the gene',
  icon: [glygen_icon],
  pagerank: 2,
})
.inputs({ protein: ProteinTerm })
.output(GlyGenProteinResponseNode)
.resolve(async (props) => {
  console.log("===> Got protein input %s", props.inputs.protein);
  const id_request = await fetch('https://api.glyGen.org/protein/search/', {
    method: 'POST',
    headers: {
      accept: 'application/json',
     'Content-Type': 'application/json',
    },
    body: JSON.stringify({ protein_name: props.inputs.protein }),
  })
  const id = await id_request.json()
  console.log("===> Got protein list %s", id);
  const protein_response = await fetch('https://api.glyGen.org/protein/list/', {
    method: 'POST',
    headers: {
      accept: 'application/json',
      'Content-Type': 'application/json',
    },
    body: JSON.stringify({ id: id['list_id'] }),
    })
  const searchResult = await protein_response.json()
  const filteredResult = filterGlyGenResults(searchResult, 'protein', props.inputs.protein);
  return filteredResult;
})
// .prompt(() => {
// })
.story(props =>
  // TODO: re-write story sentence 
  `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from.`
)
.build()

export const GlyGenProteinProduct = MetaNode('GGPP')
  .meta({
    label: 'Search GlyGen by Gene Name for Protein Products',
    description: 'Find protein product records in GlyGen for the gene',
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ gene: GeneInfo })
  .output(GlyGenProteinResponseNode)
  .resolve(async (props) => {
    const id_request = await fetch('https://api.glyGen.org/protein/search/', {
      method: 'POST',
      headers: {
        accept: 'application/json',
       'Content-Type': 'application/json',
      },
      body: JSON.stringify({ gene_name: props.inputs.gene.symbol }),
    })
    const id = await id_request.json()
    const protein_response = await fetch('https://api.glyGen.org/protein/list/', {
      method: 'POST',
      headers: {
        accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ id: id['list_id'] }),
      })
    const searchResult = await protein_response.json()
    const filteredResult = filterGlyGenResults(searchResult, 'gene', props.inputs.gene.symbol);
    return filteredResult;
  })
  // .prompt(() => {
  // })
  .story(props =>
    `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from ${props.inputs ? props.inputs.gene.symbol : 'the gene'}.`
  )
  .build()

  export const GlyGenProteinInformation = MetaNode('GlyGenProteinInformation')
  .meta({
    label: 'Search GlyGen for Protein Products',
    description: 'Find protein product records in GlyGen for the gene',
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ gene: GeneTerm })
  .output(GlyGenProteinResponseNode)
  .resolve(async (props) => {
    const gene = await GeneInfoFromGeneTerm.resolve(props)
    return await GlyGenProteinProduct.resolve({ inputs: { gene } })
  })
  .story(props =>
    `Next, Gene Info was resolved From the Gene Term "${props.inputs ? props.inputs.gene : 'the gene'}" via mygene.info.`
  )
  .build()

  // ===================================== TODO ===================================== 
  // export const GlyGenProteinQuery = MetaNode('GlyGenProteinQuery')
  // .meta({
  //   label: 'Search GlyGen by Protein Name for Protein Products',
  //   description: 'Find protein product records in GlyGen for the gene',
  //   icon: [glygen_icon],
  //   pagerank: 2,
  // })
  // .inputs({ protein_name: ProteinTerm })
  // .output(ProteinResponseNode)
  // .resolve(async (props) => {
  //   const id_request = await fetch('https://api.glygen.org/protein/search/', {
  //     method: 'POST',
  //     headers: {
  //       accept: 'application/json',
  //      'Content-Type': 'application/json',
  //     },
  //     body: JSON.stringify({ protein_name: props.inputs.protein_name }),
  //   })
  //   const id = await id_request.json()
  //   const protein_response = await fetch('https://api.glygen.org/protein/list/', {
  //     method: 'POST',
  //     headers: {
  //       accept: 'application/json',
  //       'Content-Type': 'application/json',
  //     },
  //     body: JSON.stringify({ id: id['list_id'] }),
  //     })
  //   const searchResult = await protein_response.json()
  //   const filteredResult = filterGlyGenResults(searchResult, props.inputs.protein_name);
  //   return filteredResult;
  // })
  // .prompt(() => {

  // })
  // .story(props =>
  //   `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from ${props.inputs ? props.inputs.gene.symbol : 'the gene'}.`
  // )
  // .build()

