import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { DiseaseSet, DrugSet, GeneSet, PathwaySet, PhenotypeSet, TissueSet } from '@/components/core/input/set'
import { z } from 'zod'
import { gene_icon, enrichr_icon } from '@/icons'
import { backgrounds, Disease_backgrounds, Drug_backgrounds, Pathway_backgrounds, Phenotype_backgrounds, Tissue_backgrounds, TranscriptionFactor_backgrounds } from './backgrounds'
import { DiseaseTerm, DrugTerm, GeneTerm, MetaboliteTerm, PathwayTerm, PhenotypeTerm, TissueTerm } from '@/components/core/input/term'
import * as array from '@/utils/array'
import { ScoredDiseases, ScoredDrugs, ScoredGenes, ScoredPathways, ScoredPhenotypes, ScoredTissues } from '@/components/core/input/scored'
import { Disease, Drug, Gene, Metabolite, Pathway, Phenotype, Tissue } from '@/components/core/input/primitives'

const enrichr_url = 'https://maayanlab.cloud/Enrichr'

export const EnrichrUserList = MetaNode.createData('EnrichrUserList')
  .meta({
    label: 'Enrichr User List',
    description: 'A gene set submitted to Enrichr',
    icon: [enrichr_icon, gene_icon],
  })
  .codec(z.object({
    shortId: z.string(),
    userListId: z.number(),
  }))
  .view(userlist => (
    <div className="flex-grow flex flex-row m-0" style={{ minHeight: 500 }}>
      <iframe
        className="flex-grow border-0"
        src={`${enrichr_url}/enrich?dataset=${userlist.shortId}`}
      />
    </div>
  ))
  .build()

export const EnrichrGenesetSearch = MetaNode.createProcess('EnrichrGenesetSearch')
  .meta({
    label: `Enrichr Enrichment Analysis`,
    tags: {
      'Input Type': {
        Gene: 1,
      },
      'Input Cardinality': {
        Set: 1,
      },
      Service: {
        Enrichr: 1,
      },
      'Output Type': {
        Signature: 1,
        Gene: 1,
      },
      'Output Cardinality': {
        Weighted: 1,
      },
    },
    icon: [enrichr_icon],
    description: "Fisher's exact test, odd ratio, Jaccard index",
  })
  .inputs({ geneset: GeneSet })
  .output(EnrichrUserList)
  .resolve(async (props) => {
    const formData = new FormData()
    formData.append('list', props.inputs.geneset.join('\n'))
    formData.append('description', `playbook-partnership`)
    const response = await fetch(`${enrichr_url}/addList`, {
      method: 'post',
      body: formData,
    })
    if (response.status !== 200) {
      throw new Error(`'Enrichr addList' Request failed with status ${response.status}: ${await response.text()}`)
    }
    return await response.json()
  })
  .build()

const resolveEnrichrGenesetSearchResults = async (bg: ValuesOf<typeof backgrounds>, searchResults: { shortId: string, userListId: number }) => {
  const req = await fetch(
    `${enrichr_url}/enrich?userListId=${encodeURIComponent(searchResults.userListId)}&backgroundType=${encodeURIComponent(bg.name)}`,
  )
  const res = z.object({ [bg.name]: z.array(z.tuple([
    z.number(), z.string(), z.number(), z.number(), z.number(), z.array(z.string()), z.number(), z.number(), z.number()
  ])) }).parse(await req.json())
  const results = (res[bg.name] || []).map(([
    rank, rawTerm, pvalue, zscore, combinedscore, overlapping_genes, adjusted_pvalue, unused_0, unused_1
  ]) => {
    const m = bg.termRe.exec(rawTerm)
    const term = (m && m.groups && 'term' in m.groups && m.groups.term && m.groups.term) || rawTerm
    return { term, zscore }
  })
  results.sort((a, b) => b.zscore - a.zscore)
  return results
}

export const EnrichrGenesetSearchT = [
  ...Disease_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrGenesetSearch[${bg.name}]`)
      .meta({
        label: `Extract Significant ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract significant terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrUserList })
      .output(ScoredDiseases)
      .resolve(async (props) => {
        return await resolveEnrichrGenesetSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...Drug_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrGenesetSearch[${bg.name}]`)
      .meta({
        label: `Extract Significant ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract significant terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrUserList })
      .output(ScoredDrugs)
      .resolve(async (props) => {
        return await resolveEnrichrGenesetSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...Pathway_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrGenesetSearch[${bg.name}]`)
      .meta({
        label: `Extract Significant ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract significant terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrUserList })
      .output(ScoredPathways)
      .resolve(async (props) => {
        return await resolveEnrichrGenesetSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...Phenotype_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrGenesetSearch[${bg.name}]`)
      .meta({
        label: `Extract Significant ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract significant terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrUserList })
      .output(ScoredPhenotypes)
      .resolve(async (props) => {
        return await resolveEnrichrGenesetSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...Tissue_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrGenesetSearch[${bg.name}]`)
      .meta({
        label: `Extract Significant ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract significant terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrUserList })
      .output(ScoredTissues)
      .resolve(async (props) => {
        return await resolveEnrichrGenesetSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...TranscriptionFactor_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrGenesetSearch[${bg.name}]`)
      .meta({
        label: `Extract Significant ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract significant terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrUserList })
      .output(ScoredGenes)
      .resolve(async (props) => {
        return await resolveEnrichrGenesetSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
]  

export const EnrichrGeneSearchResults = MetaNode.createData(`EnrichrGeneSearchResults`)
  .meta({
    label: `Enrichr Gene Search Results`,
    description: `Results of an Enrichr Gene Search`,
  })
  .codec(z.string())
  .view(gene => {
    return (
      <div className="flex-grow flex flex-row m-0" style={{ minHeight: 500 }}>
        <iframe
          className="flex-grow border-0"
          src={`${enrichr_url}/#find!gene=${encodeURIComponent(gene)}`}
        />
      </div>
    )
  })
  .build()

export const EnrichrGeneSearch = MetaNode.createProcess(`EnrichrGeneSearch`)
  .meta({
    label: `Enrichr Gene Search`,
    icon: [enrichr_icon],
    description: `Find terms in Enrichr Libraries containing the gene`,
  })
  .inputs({ gene: GeneTerm })
  .output(EnrichrGeneSearchResults)
  .resolve(async (props) => props.inputs.gene)
  .build()

type ValuesOf<T> = T extends Record<infer K, infer V> ? V : never

const resolveEnrichrGeneSearchResults = async (bg: ValuesOf<typeof backgrounds>, searchResults: string) => {
  const response = await fetch(
    `${enrichr_url}/genemap?json=true&gene=${encodeURIComponent(searchResults)}`,
  )
  const results = z.object({ gene: z.record(z.string(), z.array(z.string())) }).parse(await response.json())
  const terms = results.gene[bg.name] || []
  return array.unique(terms.map((term: string) => {
    const m = bg.termRe.exec(term)
    if (m && m.groups && 'term' in m.groups && m.groups.term) return m.groups.term
    else return term
  }))
}

export const EnrichrGeneSearchT = [
  ...Disease_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrGeneSearch[${bg.name}]`)
      .meta({
        label: `Extract ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrGeneSearchResults })
      .output(DiseaseSet)
      .resolve(async (props) => {
        return await resolveEnrichrGeneSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...Drug_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrGeneSearch[${bg.name}]`)
      .meta({
        label: `Extract ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrGeneSearchResults })
      .output(DrugSet)
      .resolve(async (props) => {
        return await resolveEnrichrGeneSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...Pathway_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrGeneSearch[${bg.name}]`)
      .meta({
        label: `Extract ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrGeneSearchResults })
      .output(PathwaySet)
      .resolve(async (props) => {
        return await resolveEnrichrGeneSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...Phenotype_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrGeneSearch[${bg.name}]`)
      .meta({
        label: `Extract ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrGeneSearchResults })
      .output(PhenotypeSet)
      .resolve(async (props) => {
        return await resolveEnrichrGeneSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...Tissue_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrGeneSearch[${bg.name}]`)
      .meta({
        label: `Extract ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrGeneSearchResults })
      .output(TissueSet)
      .resolve(async (props) => {
        return await resolveEnrichrGeneSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...TranscriptionFactor_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrGeneSearch[${bg.name}]`)
      .meta({
        label: `Extract ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrGeneSearchResults })
      .output(GeneSet)
      .resolve(async (props) => {
        return await resolveEnrichrGeneSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
]


export const EnrichrTermSearchResults = MetaNode.createData(`EnrichrTermSearchResults`)
  .meta({
    label: `Enrichr Term Search Results`,
    description: `Results of an Enrichr Term Search`,
  })
  .codec(z.string())
  .view(term => {
    return (
      <div className="flex-grow flex flex-row m-0" style={{ minHeight: 500 }}>
        <iframe
          className="flex-grow border-0"
          src={`${enrichr_url}/#meta!meta=${encodeURIComponent(term)}`}
        />
      </div>
    )
  })
  .build()


export const EnrichrTermTSearch = [
  { T: Disease, TermT: DiseaseTerm },
  { T: Drug, TermT: DrugTerm },
  { T: Gene, TermT: GeneTerm },
  { T: Metabolite, TermT: MetaboliteTerm },
  { T: Pathway, TermT: PathwayTerm },
  { T: Phenotype, TermT: PhenotypeTerm },
  { T: Tissue, TermT: TissueTerm },
].map(({ T, TermT }) =>
  MetaNode.createProcess(`EnrichrTermSearch[${T.name}]`)
  .meta({
    label: `Enrichr ${T.label} Term Search`,
    icon: [...T.icon, enrichr_icon],
    description: `Find ${T.label} terms in Enrichr Libraries`,
  })
  .inputs({ term: TermT })
  .output(EnrichrTermSearchResults)
  .resolve(async (props) => props.inputs.term)
  .build()
)

const resolveEnrichrTermSearchResults = async (bg: ValuesOf<typeof backgrounds>, searchResults: string) => {
  const response = await fetch(
    `${enrichr_url}/termmap?json=true&meta=${encodeURIComponent(searchResults)}`,
  )
  const results = z.object({ terms: z.record(z.string(), z.array(z.string())) }).parse(await response.json())
  const terms = results.terms[bg.name] || []
  return array.unique(terms.map((term: string) => {
    const m = bg.termRe.exec(term)
    if (m && m.groups && 'term' in m.groups && m.groups.term) return m.groups.term
    else return term
  }))
}

export const EnrichrTermSearchT = [
  ...Disease_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrTermSearch[${bg.name}]`)
      .meta({
        label: `Extract ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrTermSearchResults })
      .output(DiseaseSet)
      .resolve(async (props) => {
        return await resolveEnrichrTermSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...Drug_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrTermSearch[${bg.name}]`)
      .meta({
        label: `Extract ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrTermSearchResults })
      .output(DrugSet)
      .resolve(async (props) => {
        return await resolveEnrichrTermSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...Pathway_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrTermSearch[${bg.name}]`)
      .meta({
        label: `Extract ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrTermSearchResults })
      .output(PathwaySet)
      .resolve(async (props) => {
        return await resolveEnrichrTermSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...Phenotype_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrTermSearch[${bg.name}]`)
      .meta({
        label: `Extract ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrTermSearchResults })
      .output(PhenotypeSet)
      .resolve(async (props) => {
        return await resolveEnrichrTermSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...Tissue_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrTermSearch[${bg.name}]`)
      .meta({
        label: `Extract ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrTermSearchResults })
      .output(TissueSet)
      .resolve(async (props) => {
        return await resolveEnrichrTermSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
  ...TranscriptionFactor_backgrounds.map(bg =>
    MetaNode.createProcess(`ExtractEnrichrTermSearch[${bg.name}]`)
      .meta({
        label: `Extract ${bg.label} Signatures`,
        icon: [enrichr_icon],
        description: `Extract terms from ${bg.label} libraries`,
      })
      .inputs({ searchResults: EnrichrTermSearchResults })
      .output(GeneSet)
      .resolve(async (props) => {
        return await resolveEnrichrTermSearchResults(bg, props.inputs.searchResults)
      })
      .build()
  ),
]
