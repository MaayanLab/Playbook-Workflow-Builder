import React from 'react'
import { MetaNode, DataMetaNode, InternalDataMetaNode } from '@/spec/metanode'
import { DiseaseSet, DrugSet, GeneSet, PathwaySet, PhenotypeSet, TissueSet } from '@/components/core/input/set'
import { z } from 'zod'
import { gene_icon, enrichr_icon, search_icon } from '@/icons'
import { backgrounds, Disease_backgrounds, Drug_backgrounds, Pathway_backgrounds, Phenotype_backgrounds, Tissue_backgrounds, Gene_backgrounds } from './backgrounds'
import { DiseaseTerm, DrugTerm, GeneTerm, MetaboliteTerm, PathwayTerm, PhenotypeTerm, TissueTerm } from '@/components/core/input/term'
import { GMT } from '@/components/data/gene_matrix_transpose'
import * as array from '@/utils/array'
import * as dict from '@/utils/dict'
import { ScoredDiseases, ScoredDrugs, ScoredGenes, ScoredPathways, ScoredPhenotypes, ScoredTissues } from '@/components/core/input/scored'
import { Disease, Drug, Gene, Metabolite, Pathway, Phenotype, Tissue } from '@/components/core/primitives'
import { Table, Cell, Column } from '@/app/components/Table'
import type { ValuesOf } from '@/utils/types'
import { downloadBlob } from '@/utils/download'

const enrichr_url = 'https://maayanlab.cloud/Enrichr'

function EnrichrSet_T<T = InternalDataMetaNode>(SetT: DataMetaNode<T>) {
  return MetaNode(`Enrichr[${SetT.spec}]`)
    .meta({
      label: `Enrichr ${SetT.meta.label}`,
      description: SetT.meta.description,
      icon: [enrichr_icon, ...array.ensureArray(SetT.meta.icon)],
      color: SetT.meta.color,
    })
    .codec(z.object({ background: z.string(), terms: z.array(z.string()), set: z.array(z.string()) }))
    .view(enrichrset => (
      <Table
        height={500}
        cellRendererDependencies={[enrichrset.set]}
        numRows={enrichrset.set.length}
        enableGhostCells
        enableFocusedCell
        downloads={{
          JSON: () => downloadBlob(new Blob([JSON.stringify(enrichrset)], { type: 'application/json;charset=utf-8' }), 'data.json'),
        }}
      >
        <Column
          name={enrichrset.background}
          cellRenderer={row => <Cell key={row+''}>{enrichrset.set[row]}</Cell>}
        />
      </Table>
    ))
    .build()
}

export const EnrichrDiseaseSet = EnrichrSet_T(DiseaseSet)
export const EnrichrDrugSet = EnrichrSet_T(DrugSet)
export const EnrichrPathwaySet = EnrichrSet_T(PathwaySet)
export const EnrichrPhenotypeSet = EnrichrSet_T(PhenotypeSet)
export const EnrichrTissueSet = EnrichrSet_T(TissueSet)
export const EnrichrGeneSet = EnrichrSet_T(GeneSet)

async function resolveGenesetLibrary({ terms, background }: { background: string, terms: string[] }) {
  const req = await fetch(`https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=json&libraryName=${encodeURIComponent(background)}&term=${encodeURIComponent(terms.join(';'))}`)
  const res = await req.json()
  const gmt = z.record(z.string(), z.array(z.string())).parse(res)
  return dict.init(
    terms.map(rawTerm => {
      const m = backgrounds[background].termRe.exec(rawTerm)
      const term = (m && m.groups && 'term' in m.groups && m.groups.term && m.groups.term) || rawTerm
      return { key: rawTerm, value: { description: term, set: gmt[rawTerm] } }
    })
  )
}
export const EnrichrSetTToSetT = [
  { T: Disease, EnrichrSetT: EnrichrDiseaseSet, SetT: DiseaseSet },
  { T: Drug, EnrichrSetT: EnrichrDrugSet, SetT: DrugSet },
  { T: Pathway, EnrichrSetT: EnrichrPathwaySet, SetT: PathwaySet },
  { T: Phenotype, EnrichrSetT: EnrichrPhenotypeSet, SetT: PhenotypeSet },
  { T: Tissue, EnrichrSetT: EnrichrTissueSet, SetT: TissueSet },
  { T: Gene, EnrichrSetT: EnrichrGeneSet, SetT: GeneSet },
].flatMap(({ T, EnrichrSetT, SetT }) => [
    MetaNode(`EnrichrSetTToSetT[${T.name}]`)
      .meta({
        label: `${EnrichrSetT.meta.label} as Set`,
        icon: [enrichr_icon],
        description: `Treat Enrichr set as standard set`,
      })
      .inputs({ enrichrset: EnrichrSetT })
      .output(SetT)
      .resolve(async (props) => ({ set: props.inputs.enrichrset.set }))
      .story(props => `A ${SetT.meta.label} was extracted from the Enrichr results${props.inputs ? ` for ${props.inputs.enrichrset.background}` : ''}.`)
      .build(),
    MetaNode(`EnrichrSetTToGMT[${T.name}]`)
      .meta({
        label: `${EnrichrSetT.meta.label} as GMT`,
        icon: [enrichr_icon],
        description: `Treat Enrichr set as GMT`,
      })
      .inputs({ enrichrset: EnrichrSetT })
      .output(GMT)
      .resolve(async (props) => await resolveGenesetLibrary(props.inputs.enrichrset))
      .story(props => `A GMT was extracted from the Enrichr results${props.inputs ? ` for ${props.inputs.enrichrset.background}` : ''}.`)
      .build(),

  ]
)

export const EnrichrEnrichmentAnalysis = MetaNode('EnrichrEnrichmentAnalysis')
  .meta({
    label: 'Enrichr Enrichment Analysis',
    description: 'A gene set submitted to [Enrichr](https://maayanlab.cloud/Enrichr/)',
    icon: [enrichr_icon, gene_icon],
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
  .view(userlist => (
    <div className="flex-grow flex flex-row m-0" style={{ minHeight: 500 }}>
      {'empty' in userlist ? 
        <div className="prose">Enrichment Analysis Cannot be Performed on Empty Set</div>
        : <iframe
          className="flex-grow border-0"
          src={`${enrichr_url}/enrich?dataset=${userlist.shortId}`}
        />}
    </div>
  ))
  .build()

export const EnrichrGenesetSearch = MetaNode('EnrichrGenesetSearch')
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
    description: 'Perform Enrichment Analysis',
  })
  .inputs({ geneset: GeneSet })
  .output(EnrichrEnrichmentAnalysis)
  .resolve(async (props) => {
    if (props.inputs.geneset.set.length === 0) {
      return { empty: true }
    }
    const formData = new FormData()
    formData.append('list', props.inputs.geneset.set.join('\n'))
    formData.append('description', `playbook-partnership${props.inputs.geneset.description ? `:${props.inputs.geneset.description}` : ''}`)
    const response = await fetch(`${enrichr_url}/addList`, {
      method: 'post',
      body: formData,
    })
    if (response.status !== 200) {
      throw new Error(`'Enrichr addList' Request failed with status ${response.status}: ${await response.text()}`)
    }
    return await response.json()
  })
  .story(props =>
    `The gene set${props.inputs && props.inputs.geneset.description ? ` containing ${props.inputs.geneset.description}` : ''} was submitted to Enrichr [\\ref{doi:10.1002/cpz1.90}].`
  )
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
  { backgrounds: Disease_backgrounds, output: ScoredDiseases },
  { backgrounds: Drug_backgrounds, output: ScoredDrugs },
  { backgrounds: Pathway_backgrounds, output: ScoredPathways },
  { backgrounds: Phenotype_backgrounds, output: ScoredPhenotypes },
  { backgrounds: Tissue_backgrounds, output: ScoredTissues },
  { backgrounds: Gene_backgrounds, output: ScoredGenes },
].flatMap(({ backgrounds, output }) =>
backgrounds.map(bg => ({ bg, output }))
).map(({ bg, output }) =>
  MetaNode(`ExtractEnrichrGenesetSearch[${bg.name}]`)
    .meta({
      label: `Extract Enriched ${bg.termLabel}`,
      icon: [enrichr_icon, ...(bg.icon||[])],
      description: `Extract Significant Terms from the ${bg.label} Library`,
    })
    .inputs({ searchResults: EnrichrEnrichmentAnalysis })
    .output(output)
    .resolve(async (props) => {
      return 'empty' in props.inputs.searchResults ? [] : await resolveEnrichrGenesetSearchResults(bg, props.inputs.searchResults)
    })
    .story(props =>
      `The gene set was enriched against the ${bg.label} [${bg.ref}] library to identify statistically significant ${bg.termLabel}.`
    )
    .build()
)

export const EnrichrGeneSearchResults = MetaNode(`EnrichrGeneSearchResults`)
  .meta({
    label: `Enrichr Gene Search Results`,
    icon: [enrichr_icon, search_icon, gene_icon],
    description: `Find terms in Enrichr Libraries containing the gene`,
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

export const EnrichrGeneSearch = MetaNode(`EnrichrGeneSearch`)
  .meta({
    label: `Extract Gene Sets Containing the Gene`,
    icon: [enrichr_icon, search_icon, gene_icon],
    description: `Find terms in Enrichr Libraries containing the gene`,
  })
  .inputs({ gene: GeneTerm })
  .output(EnrichrGeneSearchResults)
  .resolve(async (props) => props.inputs.gene)
  .story(props =>
    `Gene sets containing ${props.inputs ? props.inputs.gene : 'the gene'} were queried from Enrichr [\\ref{doi:10.1002/cpz1.90}].`
  )
  .build()

const resolveEnrichrGeneSearchResults = async (bg: ValuesOf<typeof backgrounds>, searchResults: string) => {
  const response = await fetch(
    `${enrichr_url}/genemap?json=true&gene=${encodeURIComponent(searchResults)}`,
  )
  const results = z.object({ gene: z.record(z.string(), z.array(z.string())) }).parse(await response.json())
  const terms = results.gene[bg.name] || []
  return {
    background: bg.name,
    terms,
    set: array.unique(terms.map((term: string) => {
      const m = bg.termRe.exec(term)
      if (m && m.groups && 'term' in m.groups && m.groups.term) return m.groups.term
      else return term
    })),
  }
}

export const EnrichrGeneSearchT = [
  { backgrounds: Disease_backgrounds, output: EnrichrDiseaseSet },
  { backgrounds: Drug_backgrounds, output: EnrichrDrugSet },
  { backgrounds: Pathway_backgrounds, output: EnrichrPathwaySet },
  { backgrounds: Phenotype_backgrounds, output: EnrichrPhenotypeSet },
  { backgrounds: Tissue_backgrounds, output: EnrichrTissueSet },
  { backgrounds: Gene_backgrounds, output: EnrichrGeneSet },
].flatMap(({ backgrounds, output }) =>
backgrounds.map(bg => ({ bg, output }))
).map(({ bg, output }) =>
  MetaNode(`ExtractEnrichrGeneSearch[${bg.name}]`)
    .meta({
      label: `Extract ${bg.termLabel} ${bg.termAssociation} the Gene`,
      icon: [enrichr_icon, ...(bg.icon||[])],
      description: `Extract Terms from the ${bg.label} Library`,
    })
    .inputs({ searchResults: EnrichrGeneSearchResults })
    .output(output)
    .resolve(async (props) => {
      return await resolveEnrichrGeneSearchResults(bg, props.inputs.searchResults)
    })
    .story(props =>
      `Identified matching terms from the ${bg.label} [${bg.ref}] library were assembled into a collection of gene sets.`
    )
    .build()
)

export const EnrichrTermSearchResults = MetaNode(`EnrichrTermSearchResults`)
  .meta({
    label: `Enrichr Term Search Results`,
    icon: [enrichr_icon, search_icon],
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
  MetaNode(`EnrichrTermSearch[${T.name}]`)
    .meta({
      label: `Extract Gene Sets Containing the ${T.label} in the Set Label`,
      icon: [...array.ensureArray(T.icon), enrichr_icon],
      description: `Find ${T.label} Terms in Enrichr Libraries`,
    })
    .inputs({ term: TermT })
    .output(EnrichrTermSearchResults)
    .resolve(async (props) => props.inputs.term)
    .story(props =>
      `Gene sets with set labels containing ${props.inputs ? props.inputs.term : `the ${TermT.meta.label}`} were queried from Enrichr [\\ref{doi:10.1002/cpz1.90}].`
    )
    .build()
)

const resolveEnrichrTermSearchResults = async (bg: ValuesOf<typeof backgrounds>, searchResults: string) => {
  const response = await fetch(
    `${enrichr_url}/termmap?json=true&meta=${encodeURIComponent(searchResults)}`,
  )
  const results = z.object({ terms: z.record(z.string(), z.array(z.string())) }).parse(await response.json())
  const terms = results.terms[bg.name] || []
  return {
    background: bg.name,
    terms,
    set: array.unique(terms.map((term: string) => {
      const m = bg.termRe.exec(term)
      if (m && m.groups && 'term' in m.groups && m.groups.term) return m.groups.term
      else return term
    }))
  }
}

export const EnrichrTermSearchT = [
  { backgrounds: Disease_backgrounds, output: EnrichrDiseaseSet },
  { backgrounds: Drug_backgrounds, output: EnrichrDrugSet },
  { backgrounds: Pathway_backgrounds, output: EnrichrPathwaySet },
  { backgrounds: Phenotype_backgrounds, output: EnrichrPhenotypeSet },
  { backgrounds: Tissue_backgrounds, output: EnrichrTissueSet },
  { backgrounds: Gene_backgrounds, output: EnrichrGeneSet },
].flatMap(({ backgrounds, output }) =>
  backgrounds.map(bg => ({ bg, output }))
).map(({ bg, output }) =>
  MetaNode(`ExtractEnrichrTermSearch[${bg.name}]`)
    .meta({
      label: `Extract ${bg.termLabel} ${bg.termAssociation} the Term Search`,
      icon: [enrichr_icon, ...(bg.icon||[])],
      description: `Extract Terms from the ${bg.label} Library`,
      hidden: 'hidden' in bg && bg.hidden === true,
    })
    .inputs({ searchResults: EnrichrTermSearchResults })
    .output(output)
    .resolve(async (props) => {
      return await resolveEnrichrTermSearchResults(bg, props.inputs.searchResults)
    })
    .story(props =>
      `Identified matching terms from the ${bg.label} [${bg.ref}] library were assembled into a collection of gene sets.`
    )
    .build()
)
