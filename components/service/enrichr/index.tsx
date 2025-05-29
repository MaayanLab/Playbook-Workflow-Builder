import React from 'react'
import { MetaNode, DataMetaNode, InternalDataMetaNode } from '@/spec/metanode'
import {
  DiseaseSet,
  DrugSet,
  GeneSet,
  GlycanSet,
  MetaboliteSet,
  PathwaySet,
  PhenotypeSet,
  TissueSet,
} from '@/components/core/set'
import { z } from 'zod'
import { gene_icon, enrichr_icon, search_icon } from '@/icons'
import {
  backgrounds,
  Disease_backgrounds,
  Drug_backgrounds,
  Gene_backgrounds,
  Glycan_backgrounds,
  Metabolite_backgrounds,
  Pathway_backgrounds,
  Phenotype_backgrounds,
  Tissue_backgrounds,
} from './backgrounds'
import {
  DiseaseTerm,
  DrugTerm,
  GeneTerm,
  MetaboliteTerm,
  PathwayTerm,
  PhenotypeTerm,
  TissueTerm,
} from '@/components/core/term'
import { GMT } from '@/components/data/gene_matrix_transpose'
import * as array from '@/utils/array'
import * as dict from '@/utils/dict'
import {
  ScoredDiseases,
  ScoredDrugs,
  ScoredGenes,
  ScoredGlycans,
  ScoredMetabolites,
  ScoredPathways,
  ScoredPhenotypes,
  ScoredTissues,
} from '@/components/core/scored'
import { Disease, Drug, Gene, Glycan, Metabolite, Pathway, Phenotype, Tissue } from '@/components/core/primitives'
import { Table, Cell, Column } from '@/app/components/Table'
import type { ValuesOf } from '@/utils/types'
import { downloadBlob } from '@/utils/download'
import pluralize from 'pluralize'
import python from '@/utils/python'
import { PlotlyPlot } from '@/components/viz/plotly'

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
          CSV: () => downloadBlob(new Blob([
            [
              `${enrichrset.background},${SetT.meta.label}`,
              ...(enrichrset.set.map((term, i) => [JSON.stringify(enrichrset.terms[i]), JSON.stringify(term)].join(',')))
            ].join('\n')
          ], { type: 'text/csv;charset=utf-8' }), 'data.csv'),
        }}
      >
        <Column
          name={enrichrset.background}
          cellRenderer={row => <Cell key={row+''}>{enrichrset.terms[row]}</Cell>}
        />
        <Column
          name={SetT.meta.label}
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
export const EnrichrGlycanSet = EnrichrSet_T(GlycanSet)
export const EnrichrMetaboliteSet = EnrichrSet_T(MetaboliteSet)

function EnrichrScored_T<T = InternalDataMetaNode>(ScoredT: DataMetaNode<T>) {
  return MetaNode(`Enrichr[${ScoredT.spec}]`)
    .meta({
      label: `Enrichr ${ScoredT.meta.label}`,
      description: ScoredT.meta.description,
      icon: [enrichr_icon, ...array.ensureArray(ScoredT.meta.icon)],
      color: ScoredT.meta.color,
    })
    .codec(z.object({ background: z.string(), terms: z.array(z.string()), scored: z.array(z.object({ term: z.string(), zscore: z.union([z.number(), z.literal('nan'), z.literal('inf'), z.literal('-inf')]) })) }))
    .view(enrichrscored => (
      <Table
        height={500}
        cellRendererDependencies={[enrichrscored.scored]}
        numRows={enrichrscored.scored.length}
        enableGhostCells
        enableFocusedCell
        downloads={{
          JSON: () => downloadBlob(new Blob([JSON.stringify(enrichrscored)], { type: 'application/json;charset=utf-8' }), 'data.json'),
          CSV: () => downloadBlob(new Blob([
            [
              `${enrichrscored.background},${ScoredT.meta.label},ZScore`,
              ...(enrichrscored.scored.map(({ term, zscore }, i) => [JSON.stringify(enrichrscored.terms[i]), JSON.stringify(term), zscore].join(',')))
            ].join('\n')
          ], { type: 'text/csv;charset=utf-8' }), 'data.csv'),
        }}
      >
        <Column
          name={enrichrscored.background}
          cellRenderer={row => <Cell key={row+''}>{enrichrscored.terms[row]}</Cell>}
        />
        <Column
          name={ScoredT.meta.label}
          cellRenderer={row => <Cell key={row+''}>{enrichrscored.scored[row].term}</Cell>}
        />
        <Column
          name="ZScore"
          cellRenderer={row => <Cell key={row+''}>{enrichrscored.scored[row].zscore}</Cell>}
        />
      </Table>
    ))
    .build()
}

export const EnrichrScoredDiseases = EnrichrScored_T(ScoredDiseases)
export const EnrichrScoredDrugs = EnrichrScored_T(ScoredDrugs)
export const EnrichrScoredPathways = EnrichrScored_T(ScoredPathways)
export const EnrichrScoredPhenotypes = EnrichrScored_T(ScoredPhenotypes)
export const EnrichrScoredTissues = EnrichrScored_T(ScoredTissues)
export const EnrichrScoredGenes = EnrichrScored_T(ScoredGenes)
export const EnrichrScoredGlycans = EnrichrScored_T(ScoredGlycans)
export const EnrichrScoredMetabolites = EnrichrScored_T(ScoredMetabolites)

/**
 * Get the GMT, don't send too many terms per request
 */
async function resolveGmtChunked({ background, terms }: { background: string, terms: string[] }) {
  const gmt: Record<string, string[]> = {}
  for (const chunk of array.chunked(terms, 100)) {
    const req = await fetch(`https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=json&libraryName=${encodeURIComponent(background)}&term=${encodeURIComponent(chunk.join(';'))}`)
    const res = await req.json()
    Object.assign(gmt, z.record(z.string(), z.array(z.string())).parse(res))
  }
  return gmt
}

async function resolveGenesetLibrary({ terms, background }: { background: string, terms: string[] }) {
  const gmt = await resolveGmtChunked({ background, terms })
  return dict.init(
    terms.map(rawTerm => {
      const m = backgrounds[background].termRe.exec(rawTerm)
      const term = (m && m.groups && 'term' in m.groups && m.groups.term && m.groups.term) || rawTerm
      return { key: rawTerm, value: { description: term, set: gmt[rawTerm] } }
    })
  )
}

export const EnrichrTToT = [
  { T: Disease, EnrichrScoredT: EnrichrScoredDiseases, ScoredT: ScoredDiseases, EnrichrSetT: EnrichrDiseaseSet, SetT: DiseaseSet },
  { T: Drug, EnrichrScoredT: EnrichrScoredDrugs, ScoredT: ScoredDrugs, EnrichrSetT: EnrichrDrugSet, SetT: DrugSet },
  { T: Pathway, EnrichrScoredT: EnrichrScoredPathways, ScoredT: ScoredPathways, EnrichrSetT: EnrichrPathwaySet, SetT: PathwaySet },
  { T: Phenotype, EnrichrScoredT: EnrichrScoredPhenotypes, ScoredT: ScoredPhenotypes, EnrichrSetT: EnrichrPhenotypeSet, SetT: PhenotypeSet },
  { T: Tissue, EnrichrScoredT: EnrichrScoredTissues, ScoredT: ScoredTissues, EnrichrSetT: EnrichrTissueSet, SetT: TissueSet },
  { T: Gene, EnrichrScoredT: EnrichrScoredGenes, ScoredT: ScoredGenes, EnrichrSetT: EnrichrGeneSet, SetT: GeneSet },
  { T: Glycan, EnrichrScoredT: EnrichrScoredGlycans, ScoredT: ScoredGlycans, EnrichrSetT: EnrichrGlycanSet, SetT: GlycanSet },
  { T: Metabolite, EnrichrScoredT: EnrichrScoredMetabolites, ScoredT: ScoredMetabolites, EnrichrSetT: EnrichrMetaboliteSet, SetT: MetaboliteSet },
].flatMap(({ T, EnrichrScoredT, ScoredT, EnrichrSetT, SetT }) => [
  MetaNode(`EnrichrSetTToSetT[${T.name}]`)
    .meta({
      label: `${EnrichrSetT.meta.label} as Set`,
      icon: [enrichr_icon],
      description: `Load Enrichr set as standard set`,
    })
    .inputs({ enrichrset: EnrichrSetT })
    .output(SetT)
    .resolve(async (props) => ({ set: array.unique(props.inputs.enrichrset.set) }))
    .story(props => ({
      abstract: `A ${SetT.meta.label} was extracted from the Enrichr results${props.inputs?.enrichrset?.background ? ` for ${props.inputs.enrichrset.background}` : ''}.`,
      introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
      methods: `The terms in the significantly enriched gene sets from the Enrichr\\ref{doi:10.1002/cpz1.90} results in ${props.input_refs?.enrichrset} are extracted to produce ${props.output_ref}.`,
      tableLegend: `A set of terms associated with the significantly enriched gene sets from the Enrichr\\ref{doi:10.1002/cpz1.90} results in ${props.input_refs?.enrichrset}.`,
    }))
    .build(),
  MetaNode(`EnrichrSetTToGMT[${T.name}]`)
    .meta({
      label: `${EnrichrSetT.meta.label} as GMT`,
      icon: [enrichr_icon],
      description: `Load Enrichr set as GMT`,
      external: true,
    })
    .inputs({ enrichrset: EnrichrSetT })
    .output(GMT)
    .resolve(async (props) => await resolveGenesetLibrary(props.inputs.enrichrset))
    .story(props => ({
      abstract: `A GMT was extracted from the Enrichr results${props.inputs?.enrichrset?.background ? ` for ${props.inputs.enrichrset.background}` : ''}.`,
      introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
      methods: `The significantly enriched gene sets from the Enrichr\\ref{doi:10.1002/cpz1.90} results in ${props.input_refs?.enrichrset} are extracted from the original gene set library to produce ${props.output_ref}.`,
      tableLegend: `The significantly enriched gene sets filtered from the gene set library from Enrichr\\ref{doi:10.1002/cpz1.90} stored in the gene matrix transpose (GMT) format\\ref{Gene Matrix Transpose file format, https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29}.`,
    }))
    .build(),
  MetaNode(`EnrichrSetTToUMAP[${T.name}]`)
    .meta({
      label: `${EnrichrSetT.meta.label} as UMAP`,
      icon: [enrichr_icon],
      description: `Load Enrichr set as UMAP`,
      external: true,
    })
    .inputs({ enrichrset: EnrichrSetT })
    .output(PlotlyPlot)
    .resolve(async (props) => await python(
      'components.service.enrichr.resolveGenesetLibraryUMAP',
      { kargs: [props.inputs.enrichrset] },
      message => props.notify({ type: 'info', message }),
    ))
    .story(props => ({
      abstract: `Relevant terms are displayed on a gene set library UMAP\\ref{doi:10.48550/arXiv.1802.03426}${props.inputs?.enrichrset?.background ? ` for ${props.inputs.enrichrset.background}` : ''}.`,
      introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
      methods: `The significantly enriched gene sets from the Enrichr\\ref{doi:10.1002/cpz1.90} results in ${props.input_refs?.enrichrset} are highlighted in a UMAP\\ref{doi:10.48550/arXiv.1802.03426} of the original gene set library to produce ${props.output_ref}.`,
      legend: `The relevant gene sets highlighted in a UMAP\\ref{doi:10.48550/arXiv.1802.03426} of the gene set library from Enrichr\\ref{doi:10.1002/cpz1.90}.`,
    }))
    .build(),
  MetaNode(`EnrichrScoredTToScoredT[${T.name}]`)
    .meta({
      label: `${EnrichrScoredT.meta.label} as Scored`,
      icon: [enrichr_icon],
      description: `Load Enrichr scored as standard scored`,
    })
    .inputs({ enrichrscored: EnrichrScoredT })
    .output(ScoredT)
    .resolve(async (props) => props.inputs.enrichrscored.scored)
    .story(props => ({
      abstract: `A ${ScoredT.meta.label} was extracted from the Enrichr results${props.inputs?.enrichrscored?.background ? ` for ${props.inputs.enrichrscored.background}` : ''}.`,
      introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
      methods: `The terms in the significantly enriched gene sets from the Enrichr\\ref{doi:10.1002/cpz1.90} results in ${props.input_refs?.enrichrscored} are extracted to produce ${props.output_ref}.`,
      tableLegend: `A table of significantly enriched terms paired with their significance scores as reported by Enrichr\\ref{doi:10.1002/cpz1.90}.`,
    }))
    .build(),
    MetaNode(`EnrichrScoredTToUMAP[${T.name}]`)
      .meta({
        label: `${EnrichrScoredT.meta.label} as UMAP`,
        icon: [enrichr_icon],
        description: `Load Enrichr set as UMAP`,
        external: true,
      })
      .inputs({ enrichrscored: EnrichrScoredT })
      .output(PlotlyPlot)
      .resolve(async (props) => await python(
        'components.service.enrichr.resolveGeneScoredLibraryUMAP',
        { kargs: [props.inputs.enrichrscored] },
        message => props.notify({ type: 'info', message }),
      ))
      .story(props => ({
        abstract: `Relevant terms are displayed on a gene set library UMAP\\ref{doi:10.48550/arXiv.1802.03426}${props.inputs?.enrichrscored?.background ? ` for ${props.inputs.enrichrscored.background}` : ''}.`,
        introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
        methods: `The significantly enriched gene sets from the Enrichr\\ref{doi:10.1002/cpz1.90} results in ${props.input_refs?.enrichrscored} are highlighted in a UMAP\\ref{doi:10.48550/arXiv.1802.03426} of the original gene set library to produce ${props.output_ref}.`,
        legend: `The relevant gene sets highlighted in a UMAP\\ref{doi:10.48550/arXiv.1802.03426} of the gene set library from Enrichr\\ref{doi:10.1002/cpz1.90}.`,
      }))
      .build(),
  MetaNode(`EnrichrScoredTToEnrichrSetT[${T.name}]`)
    .meta({
      label: `${EnrichrScoredT.meta.label} as Enrichr Set`,
      icon: [enrichr_icon],
      description: `Load Enrichr scored as set`,
    })
    .inputs({ enrichrscored: EnrichrScoredT })
    .output(EnrichrSetT)
    .resolve(async (props) => {
      const filter = (_: unknown, i: number) => {
        const scored = props.inputs.enrichrscored.scored[i]
        return (
          scored.zscore === 'inf'
          || (typeof scored.zscore === 'number' && Math.abs(scored.zscore) > 2)
        )
      }
      return {
        background: props.inputs.enrichrscored.background,
        terms: props.inputs.enrichrscored.terms.filter(filter),
        set: props.inputs.enrichrscored.scored.filter(filter).map(({ term }) => term),
      }
    })
    .story(props => ({
      abstract: `A ${ScoredT.meta.label} was extracted from the Enrichr results${props.inputs?.enrichrscored?.background ? ` for ${props.inputs.enrichrscored.background}` : ''}.`,
      introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
      methods: `The terms in the significantly enriched gene sets from the Enrichr\\ref{doi:10.1002/cpz1.90} results in ${props.input_refs?.enrichrscored} are extracted to produce ${props.output_ref}.`,
      tableLegend: `A set of terms associated with the significantly enriched gene sets from the Enrichr\\ref{doi:10.1002/cpz1.90} results in ${props.input_refs?.enrichrscored}.`,
    }))
    .build(),
])

export const EnrichrEnrichmentAnalysis = MetaNode('EnrichrEnrichmentAnalysis')
  .meta({
    label: 'Enrichr Enrichment Analysis',
    description: 'A gene set submitted to [Enrichr](https://maayanlab.cloud/Enrichr/)',
    icon: [enrichr_icon, gene_icon],
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
  .view(userlist => (
    <div className="flex-grow flex flex-row m-0" style={{ minHeight: 500 }}>
      {'empty' in userlist ? 
        <div className="prose max-w-none">Enrichment Analysis Cannot be Performed on Empty Set</div>
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
    external: true,
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
  .story(props => ({
    abstract: `The gene set${props.inputs && props.inputs.geneset.description ? ` containing ${props.inputs.geneset.description}` : ''} was submitted to Enrichr\\ref{doi:10.1002/cpz1.90}.`,
    introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
    methods: `The gene set in ${props.input_refs?.geneset} is submitted to Enrichr\\ref{doi:10.1002/cpz1.90} for enrichment analysis.`,
    legend: `An interactive page with the results of the Enrichr\\ref{doi:10.1002/cpz1.90} enrichment analysis. Bar charts show the significantly enriched terms from different gene set libraries spanning several categories.`,
  }))
  .build()

const resolveEnrichrGenesetSearchResults = async (bg: ValuesOf<typeof backgrounds>, searchResults: { shortId: string, userListId: number }) => {
  const req = await fetch(
    `${enrichr_url}/enrich?userListId=${encodeURIComponent(searchResults.userListId)}&backgroundType=${encodeURIComponent(bg.name)}`,
  )
  const res = z.object({ [bg.name]: z.array(z.tuple([
    z.number(), z.string(), z.number(), z.number(), z.number(), z.array(z.string()), z.number(), z.number(), z.number()
  ])) }).parse(await req.json())
  const terms: string[] = []
  const scored: {term: string, zscore: number}[] = []
  for (const [
    rank, rawTerm, pvalue, zscore, combinedscore, overlapping_genes, adjusted_pvalue, unused_0, unused_1
  ] of (res[bg.name] || [])) {
    const m = bg.termRe.exec(rawTerm)
    const term = (m && m.groups && 'term' in m.groups && m.groups.term && m.groups.term) || rawTerm
    terms.push(rawTerm)
    scored.push({ term, zscore })
  }
  const sorted = array.arange(terms.length).sort((a, b) => scored[b].zscore - scored[a].zscore)
  return {
    terms: sorted.map(i => terms[i]),
    scored: sorted.map(i => scored[i]),
  }
}

export const EnrichrGenesetSearchT = [
  { backgrounds: Disease_backgrounds, output: EnrichrScoredDiseases, T: Disease },
  { backgrounds: Drug_backgrounds, output: EnrichrScoredDrugs, T: Drug },
  { backgrounds: Pathway_backgrounds, output: EnrichrScoredPathways, T: Pathway },
  { backgrounds: Phenotype_backgrounds, output: EnrichrScoredPhenotypes, T: Phenotype },
  { backgrounds: Tissue_backgrounds, output: EnrichrScoredTissues, T: Tissue },
  { backgrounds: Gene_backgrounds, output: EnrichrScoredGenes, T: Gene },
  { backgrounds: Glycan_backgrounds, output: EnrichrScoredGlycans, T: Glycan },
  { backgrounds: Metabolite_backgrounds, output: EnrichrScoredMetabolites, T: Metabolite },
].flatMap(({ backgrounds, output, T }) =>
backgrounds.map(bg => ({ bg, output, T }))
).map(({ bg, output, T }) =>
  MetaNode(`ExtractEnrichrGenesetSearch[${bg.name}]`)
    .meta({
      label: `Extract Enriched ${bg.termLabel}`,
      icon: [enrichr_icon, ...(bg.icon||[])],
      description: `Extract Significant Terms from the ${bg.label} Library`,
      tags: {
        'Output Type': {
          [pluralize(T.label)]: 1,
        },
        ...('tags' in bg ? bg.tags as Record<string, Record<string, 1>> : {}),
      },
      external: true,
    })
    .inputs({ searchResults: EnrichrEnrichmentAnalysis })
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

export const EnrichrGeneSearchResults = MetaNode(`EnrichrGeneSearchResults`)
  .meta({
    label: `Enrichr Gene Search Results`,
    icon: [enrichr_icon, search_icon, gene_icon],
    description: `Find terms in Enrichr Libraries containing the gene`,
    external: true,
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
    external: true,
  })
  .inputs({ gene: GeneTerm })
  .output(EnrichrGeneSearchResults)
  .resolve(async (props) => props.inputs.gene)
  .story(props => ({
    abstract: `Gene sets containing ${props.inputs?.gene ? props.inputs.gene : 'the gene'} were queried from Enrichr\\ref{doi:10.1002/cpz1.90}.`,
    introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
    methods: `${props.inputs?.gene ? props.inputs.gene : 'the gene'} is submitted to the Enrichr API\\ref{doi:10.1002/cpz1.90} to identify gene sets containing the gene.`,
    legend: `An interactive page provided by Enrichr\\ref{doi:10.1002/cpz1.90} showing the gene set libraries categorized by type and the gene set labels which contain the ${props.inputs?.gene ? props.inputs.gene : 'the gene'}.`,
  }))
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
    set: terms.map((term: string) => {
      const m = bg.termRe.exec(term)
      if (m && m.groups && 'term' in m.groups && m.groups.term) return m.groups.term
      else return term
    }),
  }
}

export const EnrichrGeneSearchT = [
  { backgrounds: Disease_backgrounds, output: EnrichrDiseaseSet, T: Disease },
  { backgrounds: Drug_backgrounds, output: EnrichrDrugSet, T: Drug },
  { backgrounds: Pathway_backgrounds, output: EnrichrPathwaySet, T: Pathway },
  { backgrounds: Phenotype_backgrounds, output: EnrichrPhenotypeSet, T: Phenotype },
  { backgrounds: Tissue_backgrounds, output: EnrichrTissueSet, T: Tissue },
  { backgrounds: Gene_backgrounds, output: EnrichrGeneSet, T: Gene },
  { backgrounds: Glycan_backgrounds, output: EnrichrGlycanSet, T: Glycan },
  { backgrounds: Metabolite_backgrounds, output: EnrichrMetaboliteSet, T: Metabolite },
].flatMap(({ backgrounds, output, T }) =>
backgrounds.map(bg => ({ bg, output, T }))
).map(({ bg, output, T }) =>
  MetaNode(`ExtractEnrichrGeneSearch[${bg.name}]`)
    .meta({
      label: `Extract ${bg.termLabel} ${bg.termAssociation} the Gene`,
      icon: [enrichr_icon, ...(bg.icon||[])],
      description: `Extract Terms from the ${bg.label} Library`,
      tags: {
        'Output Type': {
          [pluralize(T.label)]: 1,
        },
        ...('tags' in bg ? bg.tags as Record<string, Record<string, 1>> : {}),
      },
      external: true,
    })
    .inputs({ searchResults: EnrichrGeneSearchResults })
    .output(output)
    .resolve(async (props) => {
      return await resolveEnrichrGeneSearchResults(bg, props.inputs.searchResults)
    })
    .story(props => ({
      abstract: `Identified matching terms from the ${bg.label}${bg.ref} library were assembled into a collection of gene sets.`,
      introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
      methods: `Gene sets in the ${bg.label}${bg.ref} library containing the provided gene were resolved via the Enrichr API\\ref{doi:10.1002/cpz1.90}.`,
      tableLegend: `A table of disease signatures from the ${bg.label}${bg.ref} library from Enrichr\\ref{doi:10.1002/cpz1.90}.`,
    }))
    .build()
)

export const EnrichrTermSearchResults = MetaNode(`EnrichrTermSearchResults`)
  .meta({
    label: `Enrichr Term Search Results`,
    icon: [enrichr_icon, search_icon],
    description: `Results of an Enrichr Term Search`,
    external: true,
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
      external: true,
    })
    .inputs({ term: TermT })
    .output(EnrichrTermSearchResults)
    .resolve(async (props) => props.inputs.term)
    .story(props => ({
      abstract: `Gene sets with set labels containing ${props.inputs?.term ? props.inputs.term : `the ${TermT.meta.label}`} were queried from Enrichr\\ref{doi:10.1002/cpz1.90}.`,
      introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
      methods: `${props.inputs?.term ? props.inputs.term : `The ${TermT.meta.label}`} is submitted to the Enrichr API\\ref{doi:10.1002/cpz1.90} to identify gene sets with set labels containing the term.`,
      legend: `An interactive page provided by Enrichr\\ref{doi:10.1002/cpz1.90} showing the gene set libraries categorized by type and the gene set labels which contain ${props.inputs?.term ? props.inputs.term : 'the term'}.`,
    }))
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
    set: terms.map((term: string) => {
      const m = bg.termRe.exec(term)
      if (m && m.groups && 'term' in m.groups && m.groups.term) return m.groups.term
      else return term
    })
  }
}

export const EnrichrTermSearchT = [
  { backgrounds: Disease_backgrounds, output: EnrichrDiseaseSet, T: Disease },
  { backgrounds: Drug_backgrounds, output: EnrichrDrugSet, T: Drug },
  { backgrounds: Pathway_backgrounds, output: EnrichrPathwaySet, T: Pathway },
  { backgrounds: Phenotype_backgrounds, output: EnrichrPhenotypeSet, T: Phenotype },
  { backgrounds: Tissue_backgrounds, output: EnrichrTissueSet, T: Tissue },
  { backgrounds: Gene_backgrounds, output: EnrichrGeneSet, T: Gene },
  { backgrounds: Glycan_backgrounds, output: EnrichrGlycanSet, T: Glycan },
  { backgrounds: Metabolite_backgrounds, output: EnrichrMetaboliteSet, T: Metabolite },
].flatMap(({ backgrounds, output, T }) =>
  backgrounds.map(bg => ({ bg, output, T }))
).map(({ bg, output, T }) =>
  MetaNode(`ExtractEnrichrTermSearch[${bg.name}]`)
    .meta({
      label: `Extract ${bg.termLabel} ${bg.termAssociation} the Term Search`,
      icon: [enrichr_icon, ...(bg.icon||[])],
      description: `Extract Terms from the ${bg.label} Library`,
      hidden: 'hidden' in bg && bg.hidden === true,
      external: true,
      tags: {
        'Output Type': {
          [pluralize(T.label)]: 1,
        },
        ...('tags' in bg ? bg.tags as Record<string, Record<string, 1>> : {}),
      },
    })
    .inputs({ searchResults: EnrichrTermSearchResults })
    .output(output)
    .resolve(async (props) => {
      return await resolveEnrichrTermSearchResults(bg, props.inputs.searchResults)
    })
    .story(props => ({
      abstract: `Identified matching terms from the ${bg.label}${bg.ref} library were assembled into a collection of gene sets.`,
      introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
      methods: `Gene sets in the ${bg.label}${bg.ref} library containing the provided term in the gene set label were resolved via the Enrichr API\\ref{doi:10.1002/cpz1.90}.`,
      tableLegend: `A table of disease signatures from the ${bg.label}${bg.ref} library from Enrichr\\ref{doi:10.1002/cpz1.90}.`,
    }))
    .build()
)
