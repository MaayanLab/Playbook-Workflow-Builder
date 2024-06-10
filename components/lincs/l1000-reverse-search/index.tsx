import { z } from 'zod'
import { GeneTerm } from '@/components/core/term'
import { MetaNode } from '@/spec/metanode'
import { ScoredDrugs, ScoredGenes } from '@/components/core/scored'
import { lincs_icon, plot_icon, up_icon, down_icon } from '@/icons'

const lincs_l1000_reverse_search_dashboard = 'https://lincs-reverse-search-dashboard.dev.maayanlab.cloud'

async function lincs_l1000_reverse_search_cp({ gene, direction }: { gene: string, direction: 'up' | 'down' }) {
  const req = await fetch(`${lincs_l1000_reverse_search_dashboard}/api/table/cp/${direction}/${encodeURIComponent(gene)}`)
  const res = await req.json()
  return z.object({
    'Perturbagen': z.record(z.string(), z.string()),
    'Dose': z.record(z.string(), z.string()),
    'Timepoint': z.record(z.string(), z.string()),
    'Cell Line': z.record(z.string(), z.string()),
    'Rank in Signature': z.record(z.string(), z.number()),
    'Log2(Fold Change)': z.record(z.string(), z.number()),
    'CD Coefficient': z.record(z.string(), z.number()),
    'Fold Change': z.record(z.string(), z.number()),
  }).parse(res)
}

async function lincs_l1000_reverse_search_xpr({ gene, direction }: { gene: string, direction: 'up' | 'down' }) {
  const req = await fetch(`${lincs_l1000_reverse_search_dashboard}/api/table/xpr/${direction}/${encodeURIComponent(gene)}`)
  const res = await req.json()
  return z.object({
    'KO Gene': z.record(z.string(), z.string()),
    'Cell Line': z.record(z.string(), z.string()),
    'Timepoint': z.record(z.string(), z.string()),
    'Rank in Signature': z.record(z.string(), z.number()),
    'Log2(Fold Change)': z.record(z.string(), z.number()),
    'CD Coefficient': z.record(z.string(), z.number()),
    'Fold Change': z.record(z.string(), z.number()),
  }).parse(res)
}

export const LINCSL1000ReverseSearchDashboard = MetaNode('LINCSL1000ReverseSearchDashboard')
  .meta({
    label: `LINCS L1000 Reverse Search Dashboard`,
    description: 'A dashboard for performing L1000 Reverse Search queries for a given gene',
    icon: [lincs_icon, plot_icon],
  })
  .codec(z.object({ gene: z.string() }))
  .view(value => {
    return (
      <div className="flex-grow flex flex-row m-0" style={{ minHeight: 1100 }}>
        <iframe
          className="flex-grow border-0"
          src={`${lincs_l1000_reverse_search_dashboard}/#${encodeURIComponent(value.gene)}`}
        />
      </div>
    )
  })
  .build()

export const LINCSL1000ReverseSearch = MetaNode('LINCSL1000ReverseSearch')
  .meta({
    label: `LINCS L1000 Reverse Search`,
    description: 'Identify RNA-seq-like LINCS L1000 Signatures which reverse the expression of the gene.',
    icon: [lincs_icon],
    pagerank: 2,
  })
  .inputs({ gene: GeneTerm })
  .output(LINCSL1000ReverseSearchDashboard)
  .resolve(async (props) => {
    return { gene: props.inputs.gene }
  })
  .story(props => ({
    abstract: `RNA-seq-like LINCS L1000 Signatures [\\ref{doi:10.1093/nar/gkac328}] which mimick or reverse the the expression of ${props.inputs ? props.inputs.gene : 'the gene'} were visualized.`
  }))
  .build()

export const LINCSL1000ReverseSearchExtractDrugUp = MetaNode('LINCSL1000ReverseSearchExtract[Drug, Up]')
  .meta({
    label: `Extract Up Regulating Perturbagens`,
    description: 'Identify RNA-seq-like LINCS L1000 Chemical Perturbagen Signatures which reverse the expression of the gene.',
    icon: [up_icon],
  })
  .inputs({ search: LINCSL1000ReverseSearchDashboard })
  .output(ScoredDrugs)
  .resolve(async (props) => {
    const res = await lincs_l1000_reverse_search_cp({ gene: props.inputs.search.gene, direction: 'up' })
    const index = Object.keys(res['Perturbagen'])
    const scores = index.map(ind => res['CD Coefficient'][ind])
      .map((zscore, i) => ({ term: res['Perturbagen'][index[i]], zscore }))
    scores.sort((a, b) => b.zscore - a.zscore)
    return scores
  })
  .story(props => ({
    abstract: `Drugs which up-regulate the expression of ${props.inputs ? props.inputs.search.gene : 'the gene'} were identified from the RNA-seq-like LINCS L1000 Chemical Perturbagens [\\ref{doi:10.1093/nar/gkac328}].`
  }))
  .build()

export const LINCSL1000ReverseSearchExtractDrugDown = MetaNode('LINCSL1000ReverseSearchExtract[Drug, Down]')
  .meta({
    label: `Extract Down Regulating Perturbagens`,
    description: 'Identify RNA-seq-like LINCS L1000 Chemical Perturbagen Signatures which reverse the expression of the gene.',
    icon: [down_icon],
  })
  .inputs({ search: LINCSL1000ReverseSearchDashboard })
  .output(ScoredDrugs)
  .resolve(async (props) => {
    const res = await lincs_l1000_reverse_search_cp({ gene: props.inputs.search.gene, direction: 'down' })
    const index = Object.keys(res['Perturbagen'])
    const scores = index.map(ind => res['CD Coefficient'][ind])
      .map((zscore, i) => ({ term: res['Perturbagen'][index[i]], zscore }))
    scores.sort((a, b) => a.zscore - b.zscore)
    return scores
  })
  .story(props => ({
    abstract: `Drugs which down-regulate the expression of ${props.inputs ? props.inputs.search.gene : 'the gene'} were identified from the RNA-seq-like LINCS L1000 Chemical Perturbagens [\\ref{doi:10.1093/nar/gkac328}].`
  }))
  .build()

export const LINCSL1000ReverseSearchExtractGeneUp = MetaNode('LINCSL1000ReverseSearchExtract[Gene, Up]')
  .meta({
    label: `Extract Up Regulating CRISPR KO Genes`,
    description: 'Identify RNA-seq-like LINCS L1000 CRISPR KO Signatures which reverse the expression of the gene.',
    icon: [up_icon],
  })
  .inputs({ search: LINCSL1000ReverseSearchDashboard })
  .output(ScoredGenes)
  .resolve(async (props) => {
    const res = await lincs_l1000_reverse_search_xpr({ gene: props.inputs.search.gene, direction: 'up' })
    const index = Object.keys(res['KO Gene'])
    const scores = index.map(ind => res['CD Coefficient'][ind])
      .map((zscore, i) => ({ term: res['KO Gene'][index[i]], zscore }))
    scores.sort((a, b) => b.zscore - a.zscore)
    return scores
  })
  .story(props => ({
    abstract: `Genes which up-regulate the expression of ${props.inputs ? props.inputs.search.gene : 'the gene'} were identified from the RNA-seq-like LINCS L1000 CRISPR Knockouts [\\ref{doi:10.1093/nar/gkac328}].`
  }))
  .build()

export const LINCSL1000ReverseSearchExtractGeneDown = MetaNode('LINCSL1000ReverseSearchExtract[Gene, Down]')
  .meta({
    label: `Extract Down Regulating CRISPR KOs`,
    description: 'Identify RNA-seq-like LINCS L1000 CRISPR KO Signatures which reverse the expression of the gene.',
    icon: [down_icon],
  })
  .inputs({ search: LINCSL1000ReverseSearchDashboard })
  .output(ScoredGenes)
  .resolve(async (props) => {
    const res = await lincs_l1000_reverse_search_xpr({ gene: props.inputs.search.gene, direction: 'down' })
    const index = Object.keys(res['KO Gene'])
    const scores = index.map(ind => res['CD Coefficient'][ind])
      .map((zscore, i) => ({ term: res['KO Gene'][index[i]], zscore }))
    scores.sort((a, b) => a.zscore - b.zscore)
    return scores
  })
  .story(props => ({
    abstract: `Genes which down-regulate the expression of ${props.inputs ? props.inputs.search.gene : 'the gene'} were identified from the RNA-seq-like LINCS L1000 CRISPR Knockouts [\\ref{doi:10.1093/nar/gkac328}].`
  }))
  .build()
