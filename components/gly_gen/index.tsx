import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { GeneInfo, GeneInfoFromGeneTerm } from '../service/mygeneinfo'
import { z } from 'zod'
import { glygen_icon } from '@/icons'
import { GeneTerm, ProteinTerm, GlycanTerm } from '@/components/core/input/term'
import { filterGlyGenResults, resolveFilteredResult, GlycosylationTable, GlycanClassification, GlycanCrossRef } from './utils'

// -------- Schema Definitions -------- // 

/**
 * Zod definition for glycosylation data 
 */
const GlycosylationEntry = z.object({
  site_lbl: z.string(),
  site_category: z.string(),
  type: z.string(),
  glytoucan_ac: z.string()
}); 

/*
 * Zod definition for the overall glygen protein data  
 */
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
    species: z.object({
      name: z.string(),
      common_name: z.string(),
      taxid: z.string(),
    }),
    glycoprotein: z.object({
      glycosylation: z.boolean().optional(),
      glycosylation_data: z.array(GlycosylationEntry)
    })
  })

/**
 * Zod definition for glycan data 
 */
export const GlycanResponse = z.object({
  glytoucan: z.object({
    glytoucan_ac: z.string()
  }),
  mass: z.number(),
  mass_pme: z.number(),
  classification: z.array(
    z.object({
      type: z.object({
        name: z.string()
      }),
      subtype: z.object({
        name: z.string()
      })
    })
  ),
  crossref: z.array(
    z.object({
      id: z.string(),
      url: z.string().optional(),
      database: z.string()
    })
  )
})


// -------- Data Metanodes -------- // 

/**
 * Data metanode for the glygen api protein response, defines how the protein api response should 
 * be rendered in the UI 
 */
export const GlyGenProteinResponseNode = MetaNode('GlyGenProteinResponse')
  .meta({
    label: 'GlyGen Protein Products',
    description: 'Protein product records in GlyGen',
    icon: [glygen_icon],
  })
  .codec(GlyGenProteinResponse)
  .view( data => {
    const glyGenLink = `http://www.glygen.org/protein/${data.uniprot.uniprot_canonical_ac}`

    return (
      <div>
          <div>Gene Name: <b>{data.gene.name}</b></div>
          <div>UniProtKB Accession:
            <b>
              <a href={glyGenLink} target='_blank' rel='noopener noreferrer' style={{color: 'blue'}}> <u style={{color: 'blue'}}>{data.uniprot.uniprot_canonical_ac}</u></a>
            </b>
          </div>
          <div>Gene location: Chromosome: {data.gene.locus.chromosome} ({data.gene.locus.start_pos} - {data.gene.locus.end_pos}, '{data.gene.locus.strand}' strand)</div>
          <div>UniProtKB ID: {data.uniprot.uniprot_id}</div>
          <div>Protein Length: <b>{data.uniprot.length}</b></div>
          <div>UniProtKB Protein Name(s): {data.protein_names.name}</div>
          <div>Organism: <b>{data.species.name} ({data.species.common_name}; TaxID: {data.species.taxid})</b></div>
          <div>Glycoprotein: {data.glycoprotein.glycosylation ? 'True' : 'False'}</div>
          <br />
          <div>{data.glycoprotein.glycosylation && data.glycoprotein.glycosylation_data.length > 0 && (
            <GlycosylationTable
              glycosylationData={
                data.glycoprotein.glycosylation_data.length > 5
                ? data.glycoprotein.glycosylation_data.slice(0, 5)
                : data.glycoprotein.glycosylation_data
              }
              isPreview={data.glycoprotein.glycosylation_data.length > 5}
            />
          )}</div>
      </div>
    )
  })
  .build()

/**
 * Data metanode for the glycosylation list from a protein response, defines how the glycosylation 
 * information contained within the protein response json should be rendered in the UI
 */
export const GlycosylationViewResponseNode = MetaNode('GlycosylationViewResponse')
  .meta({
    label: 'Glycosylation Information for Glycoproteins',
    description: 'Glycosylation product records in GlyGen',
    icon: [glygen_icon],
  })
  .codec(GlyGenProteinResponse)
  .view( data => {
    const recordCount = data.glycoprotein.glycosylation_data.length
    const emptyGlytoucanAcCount = data.glycoprotein.glycosylation_data.filter(entry => entry.glytoucan_ac === '').length
    const glyGenLink = `http://www.glygen.org/protein/${data.uniprot.uniprot_canonical_ac}`

    if (recordCount > 0) {
      return (
        <div>
          <div>
            UniProtKB Accession: 
              <b>
                <a href={glyGenLink} target='_blank' rel='noopener noreferrer' style={{color: 'blue'}}> <u style={{color: 'blue'}}>{data.uniprot.uniprot_canonical_ac}</u></a>
              </b>
          </div>
          <div>Total Records: <b>{data.glycoprotein.glycosylation_data.length}</b></div>
          <div>Records Without GlyTouCan Accessions: <b>{emptyGlytoucanAcCount}</b></div>
          <br/>
          <GlycosylationTable
            glycosylationData={data.glycoprotein.glycosylation_data}
          />
        </div>
      )
    } else {
      return (
        <div>
          <div>
            UniProtKB Accession: 
                <b>
                  <a href={glyGenLink} target='_blank' rel='noopener noreferrer' style={{color: 'blue'}}> <u style={{color: 'blue'}}>{data.uniprot.uniprot_canonical_ac}</u></a>
                </b>
          </div>
          <div><b>No Glycosylation Information to Display</b></div>
        </div>
      )
    }
  })
  .build()

export const GlycanViewResponseNode = MetaNode('GlycanViewResponse')
  .meta({
    label: 'Glycan information',
    description: 'Glycan information from GlyGen'
  })
  .codec(GlycanResponse)
  .view(data => {
    const glyGenLink = `http://www.glygen.org/glycan/${data.glytoucan.glytoucan_ac}`

    return (
      <div>
        <div>GlyTouCan Accession: 
          <b>
            <a href={glyGenLink} target='_blank' rel='noopener nonreferrer' style={{color: 'blue'}}> <u style={{color: 'blue'}}>{data.glytoucan.glytoucan_ac}</u></a>
          </b>
        </div>
        <div>Monoisotopic Mass: <b>{data.mass} Da</b></div>
        <div>Monoisotopic Mass-pMe (Da): <b>{data.mass_pme} Da</b></div>
        <div>
          Glycan Type / Glycan Subtype: <b><GlycanClassification classification={data.classification}/></b>
        </div>
        <div>
          <GlycanCrossRef crossref={data.crossref}/>
        </div>
        <div>
          Glycan Image: 
          <img src={`https://api.glygen.org/glycan/image/${data.glytoucan.glytoucan_ac}/`} alt='Glycan Image'/>
        </div>
      </div>
    )
  })
  .build()


// -------- Process Metanodes (all resolver process metanodes) -------- // 

/**
 * Process metanode for searching by protein name for protein products 
 */
export const GlyGenProtein = MetaNode('GGP')
  .meta({
    label: 'Search GlyGen by Protein Name for Protein Products',
    description: 'Find protein product records in GlyGen for the gene',
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ protein_uniprot_canonical_ac: ProteinTerm })
  .output(GlyGenProteinResponseNode)
  .resolve(async (props) => {
    const protein_response = await resolveFilteredResult(props.inputs.protein_uniprot_canonical_ac);
    return protein_response;
  })
  .story(props =>
    // TODO: re-write story sentence to make sense with protein term input (previous gene value removed to prevent `npm run build` error)
    `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from.`
  )
  .build()

/**
 * Process metanode for searching the glygen database by gene name for protein products 
 * given a GeneInfo 
 */
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
  .story(props =>
    `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from ${props.inputs ? props.inputs.gene.symbol : 'the gene'}.`
  )
  .build()

/**
 * Process metanode for searching the glygen database by gene name for protein products 
 * given a GeneTerm 
 */
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
    return await GlyGenProteinProduct.resolve({ ...props, inputs: { gene } })
  })
  .story(props =>
    `The GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of protein products that originate from ${props.inputs ? props.inputs.gene : 'the gene'}.`
  )
  .build()

/**
 * Process metanode to extract the glycosylation information from the protein prodct 
 */
export const GlycosylationInformation = MetaNode('GlycosylationInformation')
  .meta({
    label: 'Get Glycosylation Data from GlyGen Protein Products',
    description: 'Glycosylation information for Glycoproteins',
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ glygenProteinResponse: GlyGenProteinResponseNode })
  .output(GlycosylationViewResponseNode)
  .resolve(async (props) => {
    return props.inputs.glygenProteinResponse;
  })
  .story( props => 
    'The glycosylation data was extracted from the GlyGen protein response and prepared for presentation in the view metanode.'
  )
  .build() 

/**
 * Process metanode to search GlyGen for for Glycan details given a Glytoucan Accession
 */
export const GlycanInformation = MetaNode('GlycanInformation')
  .meta({
    label: 'Search GlyGen by GlyTouCan Accession',
    description: 'Search for Glycan information',
    // icon: []
    pagerank: 2
  })
  .inputs({ glycan: GlycanTerm })
  .output(GlycanViewResponseNode)
  .resolve(async (props) => {
    // get glycan data 
    const detail_response = await fetch(`https://api.glygen.org/glycan/detail/${props.inputs.glycan}`, {
      method: 'POST',
      headers: {
        accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ glytoucan_ac: props.inputs.glycan })
    })
    const glycan_data = await detail_response.json();
    return glycan_data
  })
  .story(props => 
    `The GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a information about ${props.inputs ? props.inputs.glycan : 'the glycan'}.`
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

// export type GlyGenProteinResponseType = z.infer<typeof GlyGenProteinResponse>

 // mass: z.object({
    //   chemical_mass: z.string()
    // }),
    // refseq_ac: z.string(),
    // refseq_name: z.string(),

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