import { MetaNode } from '@/spec/metanode'
import { ScoredGenes } from '@/components/core/scored'
import { archs4_icon, enrichr_icon, geneshot_icon } from '@/icons'
import { GeneSet } from '@/components/core/set'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import * as math from '@/utils/math'
import { Disease, Drug, Gene, Pathway, Phenotype, Tissue } from '@/components/core/primitives'
import { DiseaseTerm, DrugTerm, GeneTerm, PathwayTerm, PhenotypeTerm, TissueTerm } from '@/components/core/term'

export const geneshot_url = 'https://maayanlab.cloud/geneshot'

async function geneshot_term_search(body: { rif: 'generif' | 'autorif', term: string }) {
  const req = await fetch(`${geneshot_url}/api/search`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Accept': 'application/json',
    },
    body: JSON.stringify(body),
  })
  if (!req.ok) throw new Error(`Failed to query Geneshot: ${await req.text()}`)
  const res = z.object({
    PubMedID_count: z.number(),
    gene_count: z.record(z.string(), z.tuple([
      z.number().describe('Total number of matching PubMed ids'),
      z.number().describe('Fraction of matching PubMed ids'),
    ])),
    query_time: z.number(),
    return_size: z.number(),
    search_term: z.string(),
  }).parse(await req.json())
  const { mu, std } = math.mean_std(dict.values(res.gene_count).map(([total, fraction]) => total))
  const output = dict.items(res.gene_count).map(({ key: term, value: [total, fraction] }) => ({
    term,
    zscore: (total - mu)/std
  }))
  output.sort((a, b) => b.zscore - a.zscore)
  return output
}

export const GeneshotTermSearchT = [
  { rc: 'generif', label: 'GeneRIF', icon: [] } as const,
  { rc: 'autorif', label: 'AutoRIF', icon: [] } as const,
].flatMap(({ rc, label, icon }) =>
  [
    { T: Disease, TermT: DiseaseTerm },
    { T: Drug, TermT: DrugTerm },
    { T: Gene, TermT: GeneTerm },
    { T: Pathway, TermT: PathwayTerm },
    { T: Phenotype, TermT: PhenotypeTerm },
    { T: Tissue, TermT: TissueTerm },
  ].flatMap(({ T, TermT }) => [
    MetaNode(`GeneshotTermSearch[${T.name}, ${rc}]`)
      .meta({
        label: `Search PubMed for Gene Co-Mentions with ${label}`,
        description: `Geneshot ${T.label.toLowerCase()}-gene co-mentions`,
        icon: [geneshot_icon, ...icon],
        tags: {
          'Input Type': {
            [T.name]: 1,
          },
          'Input Cardinality': {
            Term: 1,
          },
          Service: {
            Geneshot: 1,
          },
          'Output Type': {
            Gene: 1,
          },
          'Output Cardinality': {
            Weighted: 1,
          },
        },
      })
      .inputs({ term: TermT })
      .output(ScoredGenes)
      .resolve(async (props) => {
        return await geneshot_term_search({ rif: rc, term: props.inputs.term })
      })
      .story(props => ({
        abstract: `${props.inputs ? props.inputs.term : `The ${T.label.toLowerCase()}`}-gene co-mentions on PubMed were queried with Geneshot based on ${label} [\\ref{doi:10.1093/nar/gkz393}].`
      }))
      .build()
  ])
)

async function geneshot_geneset_augmentation(body: { gene_list: string[], similarity: 'generif' | 'tagger' | 'coexpression' | 'enrichr' }) {
  const req = await fetch(`${geneshot_url}/api/associate`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Accept': 'application/json',
    },
    body: JSON.stringify(body),
  })
  if (!req.ok) throw new Error(`Failed to augment gene set using Geneshot: ${await req.text()}`)
  const res = z.object({
    association: z.record(z.string(), z.object({
      publications: z.number(),
      simScore: z.number(),
    })),
    darkgpcr: z.array(z.string()),
    darkionchannel: z.array(z.string()),
    darkkinase: z.array(z.string()),
    gpcr: z.array(z.string()),
    ionchannel: z.array(z.string()),
    kinase: z.array(z.string()),
  }).parse(await req.json())
  const output = dict.items(res.association).map(({ key: term, value }) => ({
    term,
    zscore: Number(value.simScore),
  }))
  output.sort((a, b) => b.zscore - a.zscore)
  return output
}

export const GeneshotGeneSetAugmentation = ([
  { rc: 'generif', label: 'Literature Co-Mentions with GeneRIF', icon: [] } as const,
  { rc: 'tagger', label: 'Literature Co-Mentions with Tagger', icon: [] } as const,
  { rc: 'coexpression', label: 'mRNA Co-Expression from ARCHS4', icon: [archs4_icon] } as const,
  { rc: 'enrichr', label: 'Enrichr Query Co-Occurrence', icon: [enrichr_icon] } as const,
]).flatMap(({ rc, label, icon }) => [
  MetaNode(`GeneshotGeneSetAugmentation[${rc}]`)
    .meta({
      label: `Expand a Gene Set based on ${label}`,
      description: 'Geneshot Gene Set Augmentation',
      icon: [geneshot_icon, ...icon],
      tags: {
        'Input Type': {
          Gene: 1,
        },
        'Input Cardinality': {
          Set: 1,
        },
        Service: {
          Geneshot: 1,
        },
        'Output Type': {
          Gene: 1,
        },
        'Output Cardinality': {
          Weighted: 1,
        },
      },
    })
    .inputs({ geneset: GeneSet })
    .output(ScoredGenes)
    .resolve(async (props) => {
      return await geneshot_geneset_augmentation({
        similarity: rc,
        gene_list: props.inputs.geneset.set,
      })
    })
    .story(props => ({
      abstract: `The gene set${props.inputs && props.inputs.geneset.description ? ` containing ${props.inputs.geneset.description}` : ''} was augmented with Geneshot based on ${label} [\\ref{doi:10.1093/nar/gkz393}].`
    }))
    .build(),
  MetaNode(`GeneshotGeneAugmentation[${rc}]`)
    .meta({
      label: `Construct a Gene Set based on ${label}`,
      description: 'Geneshot Gene Set Augmentation',
      icon: [geneshot_icon, ...icon],
      tags: {
        'Input Type': {
          Gene: 1,
        },
        'Input Cardinality': {
          Term: 1,
        },
        Service: {
          Geneshot: 1,
        },
        'Output Type': {
          Gene: 1,
        },
        'Output Cardinality': {
          Weighted: 1,
        },
      },
    })
    .inputs({ gene: GeneTerm })
    .output(ScoredGenes)
    .resolve(async (props) => {
      return await geneshot_geneset_augmentation({
        similarity: rc,
        gene_list: [props.inputs.gene],
      })
    })
    .story(props => ({
      abstract: `${props.inputs ? props.inputs.gene : 'The gene'} was augmented with Geneshot based on ${label} [\\ref{doi:10.1093/nar/gkz393}].`
    }))
    .build(),
])
