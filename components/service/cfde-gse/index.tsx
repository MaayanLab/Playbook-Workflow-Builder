import { MetaNode } from "@/spec/metanode";
import { drc_icon } from "@/icons"
import { z } from 'zod'
import pluralize from "pluralize"
import * as str from '@/utils/str'
import { GeneSet } from "@/components/core/set"
import { Disease, Drug, Gene, Glycan, Metabolite, Phenotype, Tissue } from "@/components/core/primitives"
import { Disease_backgrounds, Drug_backgrounds, Phenotype_backgrounds, Tissue_backgrounds, Gene_backgrounds, Glycan_backgrounds, Metabolite_backgrounds } from "@/components/service/enrichr/backgrounds"
import { EnrichrEnrichmentAnalysis, EnrichrScoredDiseases, EnrichrScoredDrugs, EnrichrScoredGenes, EnrichrScoredGlycans, EnrichrScoredMetabolites, EnrichrScoredPathways, EnrichrScoredPhenotypes, EnrichrScoredTissues, resolveEnrichrGenesetSearchResults } from "@/components/service/enrichr"

const cfde_gse_url = 'https://gse.cfde.cloud'
const cfde_gse_enrich_url = 'https://maayanlab.cloud/Enrichr'
const cfde_gse_icon = drc_icon

export const CFDEGSEKG = MetaNode('CFDEGSEKG')
  .meta({
    label: 'CFDE Gene Set Enrichment Knowledge Graph',
    description: 'A CFDE-Centric knowledge graph for gene set enrichment results',
    icon: [cfde_gse_icon],
    external: true,
  })
  .codec(z.union([
    z.object({
      empty: z.literal(true),
    }),
    z.object({
      shortId: z.string(),
      userListId: z.number(),
    })
  ]))
  .view(userlist => {
    const searchParams = new URLSearchParams()
    if (!('empty' in userlist)) {
      searchParams.append('q', JSON.stringify({
        "min_lib": 1,
        "libraries":[
          {"name":"LINCS_L1000_Chem_Pert_Consensus_Sigs","limit":5},
          {"name":"HuBMAP_ASCTplusB_augmented_2022","limit":5},
        ],
        "userListId":userlist.userListId.toString(),
        "search":true,
      }))
      searchParams.append('fullscreen', 'true')
    }
    return (
      <div className="flex-grow flex flex-row m-0" style={{ minHeight: 1000 }}>
        {'empty' in userlist ? 
          <div className="prose max-w-none">Enrichment Analysis Cannot be Performed on Empty Set</div>
          : <iframe
            className="flex-grow border-0"
            src={`${cfde_gse_url}?${searchParams.toString()}`}
          />}
      </div>
    )
  })
  .build()

export const EnrichrUserListToCFDEGSEKG = MetaNode('EnrichrUserListToCFDEGSEKG')
  .meta({
    label: 'Visualize CFDE Gene Set Enrichment Knowledge Graph',
    description: 'View gene set results as a knowledge graph',
    tags: {
      'Output Type': {
        'Knowledge Graph': 1,
      },
      'Data Source': {
        'CFDE': 1
      },
    },
    external: true,
  })
  .inputs({ enrichr: EnrichrEnrichmentAnalysis })
  .output(CFDEGSEKG)
  .resolve(async (props) => props.inputs.enrichr)
  .story(props => ({
    introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
    methods: `The enrichment results from Enrichr\\ref{doi:10.1002/cpz1.90} are used to subset a Knowledge Graph of inter-connected CFDE Gene Set Libraries using CFDE GSE\\ref{CFDE GSE, https://gse.cfde.cloud/}.`,
    legend: `An interactive page containing a knowledge graph highlighting the enriched terms in the CFDE Gene Set Libraries provided by CFDE GSE\\ref{CFDE GSE, https://gse.cfde.cloud/}.`,
  }))
  .build()

export const CFDEGSEGenesetSearch = MetaNode('CFDEGSEGenesetSearch')
  .meta({
    label: `CFDE Enrichment Analysis`,
    tags: {
      'Input Type': {
        Gene: 1,
      },
      'Input Cardinality': {
        Set: 1,
      },
      Service: {
        'CFDE GSE': 1,
      },
      'Output Type': {
        Signature: 1,
        Gene: 1,
      },
      'Output Cardinality': {
        Weighted: 1,
      },
    },
    icon: [cfde_gse_icon],
    description: 'Perform Enrichment Analysis',
    external: true,
  })
  .inputs({ geneset: GeneSet })
  .output(CFDEGSEKG)
  .resolve(async (props) => {
    if (props.inputs.geneset.set.length === 0) {
      return { empty: true }
    }
    const formData = new FormData()
    formData.append('list', props.inputs.geneset.set.join('\n'))
    formData.append('description', str.truncate(`playbook-partnership${props.inputs.geneset.description ? `:${props.inputs.geneset.description}` : ''}`, 255))
    const response = await fetch(`${cfde_gse_enrich_url}/addList`, {
      method: 'post',
      body: formData,
    })
    if (response.status !== 200) {
      throw new Error(`'CFDE GSE addList' Request failed with status ${response.status}: ${await response.text()}`)
    }
    return await response.json()
  })
  .story(props => ({
    abstract: `The gene set${props.inputs && props.inputs.geneset.description ? ` containing ${props.inputs.geneset.description}` : ''} was submitted to CFDE GSE\\ref{Common Fund Data Ecosystem (CFDE) Gene Set Enrichment (GSE), https://gse.cfde.cloud}.`,
    introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. CFDE GSE is a gene set search engine & knowledge graph that enables the querying of CF annotated gene sets\\ref{Common Fund Data Ecosystem (CFDE) Gene Set Enrichment (GSE), https://gse.cfde.cloud}.`,
    methods: `The gene set in ${props.input_refs?.geneset} is submitted to CFDE GSE\\ref{Common Fund Data Ecosystem (CFDE) Gene Set Enrichment (GSE), https://gse.cfde.cloud} for enrichment analysis.`,
    legend: `An interactive page with the results of the CFDE GSE\\ref{Common Fund Data Ecosystem (CFDE) Gene Set Enrichment (GSE), https://gse.cfde.cloud} enrichment analysis results.`,
  }))
  .build()

export const CFDEGSEGenesetSearchT = [
  { backgrounds: Disease_backgrounds.filter(bg => 'tags' in bg && bg.tags && 'Data Source' in bg.tags && bg.tags['Data Source']['CFDE'] === 1),
    output: EnrichrScoredDiseases, T: Disease },
  { backgrounds: Drug_backgrounds.filter(bg => 'tags' in bg && bg.tags && 'Data Source' in bg.tags && bg.tags['Data Source']['CFDE'] === 1),
    output: EnrichrScoredDrugs, T: Drug },
  { backgrounds: Phenotype_backgrounds.filter(bg => 'tags' in bg && bg.tags && 'Data Source' in bg.tags && bg.tags['Data Source']['CFDE'] === 1),
    output: EnrichrScoredPhenotypes, T: Phenotype },
  { backgrounds: Tissue_backgrounds.filter(bg => 'tags' in bg && bg.tags && 'Data Source' in bg.tags && bg.tags['Data Source']['CFDE'] === 1),
    output: EnrichrScoredTissues, T: Tissue },
  { backgrounds: Gene_backgrounds.filter(bg => 'tags' in bg && bg.tags && 'Data Source' in bg.tags && bg.tags['Data Source']['CFDE'] === 1),
    output: EnrichrScoredGenes, T: Gene },
  { backgrounds: Glycan_backgrounds.filter(bg => 'tags' in bg && bg.tags && 'Data Source' in bg.tags && bg.tags['Data Source']['CFDE'] === 1),
    output: EnrichrScoredGlycans, T: Glycan },
  { backgrounds: Metabolite_backgrounds.filter(bg => 'tags' in bg && bg.tags && 'Data Source' in bg.tags && bg.tags['Data Source']['CFDE'] === 1),
    output: EnrichrScoredMetabolites, T: Metabolite },
].flatMap(({ backgrounds, output, T }) =>
backgrounds.map(bg => ({ bg, output, T }))
).map(({ bg, output, T }) =>
  MetaNode(`ExtractCFDEGSEGenesetSearch[${bg.name}]`)
    .meta({
      label: `Extract Enriched ${bg.termLabel}`,
      icon: [cfde_gse_icon, ...(bg.icon||[])],
      description: `Extract Significant Terms from the ${bg.label} Library`,
      tags: {
        'Output Type': {
          [pluralize(T.label)]: 1,
        },
        ...('tags' in bg ? bg.tags as Record<string, Record<string, 1>> : {}),
      },
      external: true,
    })
    .inputs({ searchResults: CFDEGSEKG })
    .output(output)
    .resolve(async (props) => {
      return !('empty' in props.inputs.searchResults) ? {
        background: bg.name,
        ...(await resolveEnrichrGenesetSearchResults(bg, props.inputs.searchResults))
      } : {
        background: bg.name,
        terms: [],
        scored: [],
      }
    })
    .story(props => ({
      abstract: `The gene set was enriched against the ${bg.label}${bg.ref} library to identify statistically significant ${bg.termLabel}.`,
      introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
      methods: `The enrichment analysis results against the ${bg.label}${bg.ref} library are resolved via the Enrichr API\\ref{doi:10.1002/cpz1.90}.`,
      tableLegend: `A table of significantly enriched ${bg.termLabel} from the ${bg.label}${bg.ref} library paired with their significance scores as reported by Enrichr\\ref{doi:10.1002/cpz1.90}.`,
    }))
    .build()
)
