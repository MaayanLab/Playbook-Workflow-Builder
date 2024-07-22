import { MetaNode } from '@/spec/metanode'
import { ScoredDrugs } from '@/components/core/scored'
import { archs4_icon, drugshot_icon } from '@/icons'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import * as math from '@/utils/math'
import { Disease, Drug, Gene, Pathway, Phenotype, Tissue } from '@/components/core/primitives'
import { DiseaseTerm, DrugTerm, GeneTerm, PathwayTerm, PhenotypeTerm, TissueTerm } from '@/components/core/term'
import { DrugSet } from '@/components/core/set'
import Citable from '@/utils/citations'

export const drugshot_url = 'https://maayanlab.cloud/drugshot'

async function drugshot_term_search(body: { rif: 'drugrif' | 'autorif', term: string }) {
  const req = await fetch(`${drugshot_url}/api/search`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Accept': 'application/json',
    },
    body: JSON.stringify(body),
  })
  if (!req.ok) throw new Error(`Failed to query DrugShot: ${await req.text()}`)
  const res = z.object({
    PubMedID_count: z.number(),
    drug_count: z.record(z.string(), z.tuple([
      z.number().describe('Total number of matching PubMed ids'),
      z.number().describe('Fraction of matching PubMed ids'),
    ])),
    query_time: z.number(),
    return_size: z.number(),
    search_term: z.string(),
  }).parse(await req.json())
  const { mu, std } = math.mean_std(dict.values(res.drug_count).map(([total, fraction]) => total))
  const output = dict.items(res.drug_count).map(({ key: term, value: [total, fraction] }) => ({
    term,
    zscore: (total - mu)/std
  }))
  output.sort((a, b) => b.zscore - a.zscore)
  return output
}

export const DrugShotTermSearchT = [
  { rc: 'drugrif', label: 'DrugRIF', icon: [] } as const,
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
    MetaNode(`DrugShotTermSearch[${T.name}, ${rc}]`)
      .meta({
        label: `Search PubMed for Drug Co-Mentions with ${label}`,
        description: `DrugShot ${T.label.toLowerCase()}-drug co-mentions`,
        icon: [drugshot_icon, ...icon],
        tags: {
          'Input Type': {
            [T.name]: 1,
          },
          'Input Cardinality': {
            Term: 1,
          },
          Service: {
            DrugShot: 1,
          },
          'Output Type': {
            Drug: 1,
          },
          'Output Cardinality': {
            Weighted: 1,
          },
        },
      })
      .inputs({ term: TermT })
      .output(ScoredDrugs)
      .resolve(async (props) => {
        return await drugshot_term_search({ rif: rc, term: props.inputs.term })
      })
      .story(props => ({
        abstract: Citable.text`${props.inputs?.term ? props.inputs.term : `The ${T.label.toLowerCase()}`}-drug co-mentions on PubMed were queried with DrugShot based on ${label} [${Citable.doi('10.1186/s12859-022-04590-5')}].`
      }))
      .build()
  ])
)

async function drugshot_drugset_augmentation(body: { drug_list: string[], similarity: 'drugrif_cooccur' | 'L1000_coexpression' }) {
  const req = await fetch(`${drugshot_url}/api/associate`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Accept': 'application/json',
    },
    body: JSON.stringify(body),
  })
  if (!req.ok) throw new Error(`Failed to augment drug set using DrugShot: ${await req.text()}`)
  const res = z.object({
    association: z.record(z.string(), z.object({
      publications: z.number(),
      simScore: z.number(),
    })),
  }).parse(await req.json())
  const output = dict.items(res.association).map(({ key: term, value }) => ({
    term,
    zscore: Number(value.simScore),
  }))
  output.sort((a, b) => b.zscore - a.zscore)
  return output
}

export const DrugShotDrugSetAugmentation = ([
  { rc: 'drugrif_cooccur', label: 'Literature Co-Mentions with DrugRIF', icon: [] } as const,
  { rc: 'L1000_coexpression', label: 'Co-Expression from LINCS L1000', icon: [archs4_icon] } as const,
]).flatMap(({ rc, label, icon }) => [
  MetaNode(`DrugShotDrugSetAugmentation[${rc}]`)
    .meta({
      label: `Expand a Drug Set based on ${label}`,
      description: 'DrugShot Drug Set Augmentation',
      icon: [drugshot_icon, ...icon],
      tags: {
        'Input Type': {
          Drug: 1,
        },
        'Input Cardinality': {
          Set: 1,
        },
        Service: {
          DrugShot: 1,
        },
        'Output Type': {
          Drug: 1,
        },
        'Output Cardinality': {
          Weighted: 1,
        },
      },
    })
    .inputs({ drugset: DrugSet })
    .output(ScoredDrugs)
    .resolve(async (props) => {
      return await drugshot_drugset_augmentation({
        similarity: rc,
        drug_list: props.inputs.drugset.set,
      })
    })
    .story(props => ({
      abstract: Citable.text`The drug set${props.inputs && props.inputs.drugset.description ? ` containing ${props.inputs.drugset.description}` : ''} was augmented with DrugShot based on ${label} [${Citable.doi('10.1186/s12859-022-04590-5')}].`
    }))
    .build(),
  MetaNode(`DrugShotDrugAugmentation[${rc}]`)
    .meta({
      label: `Construct a Drug Set based on ${label}`,
      description: 'DrugShot Drug Set Augmentation',
      icon: [drugshot_icon, ...icon],
      tags: {
        'Input Type': {
          Drug: 1,
        },
        'Input Cardinality': {
          Term: 1,
        },
        Service: {
          DrugShot: 1,
        },
        'Output Type': {
          Drug: 1,
        },
        'Output Cardinality': {
          Weighted: 1,
        },
      },
    })
    .inputs({ drug: DrugTerm })
    .output(ScoredDrugs)
    .resolve(async (props) => {
      return await drugshot_drugset_augmentation({
        similarity: rc,
        drug_list: [props.inputs.drug],
      })
    })
    .story(props => ({
      abstract: Citable.text`${props.inputs?.drug ? props.inputs.drug : 'The drug'} was augmented with DrugShot based on ${label} [${Citable.doi('10.1186/s12859-022-04590-5')}].`
    }))
    .build(),
])
