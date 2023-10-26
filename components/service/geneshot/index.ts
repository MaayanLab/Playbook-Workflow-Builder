import { MetaNode } from '@/spec/metanode'
import { ScoredGenes } from '@/components/core/input/scored'
import { archs4_icon, enrichr_icon, geneshot_icon } from '@/icons'
import { GeneSet } from '@/components/core/input/set'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import { GeneTerm } from '@/components/core/input/term'

export const geneshot_url = 'https://maayanlab.cloud/geneshot'

async function geneshot_geneset_augmentation(body: { gene_list: string[], similarity: string }) {
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

export const GeneshotGeneSetAugmentation = [
  { rc: 'generif', label: 'Literature Co-Mentions with GeneRIF', icon: [] },
  { rc: 'tagger', label: 'Literature Co-Mentions with Tagger', icon: [] },
  { rc: 'coexpression', label: 'mRNA Co-Expression from ARCHS4', icon: [archs4_icon] },
  { rc: 'enrichr', label: 'Enrichr Query Co-Occurrence', icon: [enrichr_icon] },
].flatMap(({ rc, label, icon }) => [
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
    .story(props =>
      `The gene set${props.inputs && props.inputs.geneset.description ? ` containing ${props.inputs.geneset.description}` : ''} was augmented with Geneshot based on ${label} [\\ref{doi:10.1093/nar/gkz393}].`
    )
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
    .story(props =>
      `${props.inputs ? props.inputs.gene : 'The gene'} was augmented with Geneshot based on ${label} [\\ref{doi:10.1093/nar/gkz393}].`
    )
    .build(),
])
