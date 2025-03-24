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
    external: true,
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
    external: true,
  })
  .inputs({ gene: GeneTerm })
  .output(LINCSL1000ReverseSearchDashboard)
  .resolve(async (props) => {
    return { gene: props.inputs.gene }
  })
  .story(props => ({
    abstract: `RNA-seq-like LINCS L1000 Signatures\\ref{doi:10.1093/nar/gkac328} which maximally increase or decrease the expression of ${props.inputs?.gene ? props.inputs.gene : 'the gene'} were visualized.`,
    introduction: `The LINCS L1000 data\\ref{doi:10.1016/j.cell.2017.10.049} now consists of over 3 million chemical and genetic perturbational signatures. A modified CycleGAN model was used to transorm the L1000 assay measurements to RNA-seq-like expression\\ref{doi:10.1186/s12859-022-04895-5}. Characteristic direction\\ref{doi:10.1186/1471-2105-15-79} was used to compute differential gene expression signatures.`,
    methods: `Pre-computed characteristic direction\\ref{doi:10.1186/1471-2105-15-79}-based differential gene expression signatures from RNA-seq-like LINCS L1000 Expression Profiles\\ref{doi:10.1186/s12859-022-04895-5} which maximally increase or decrease the expression of the gene are retrieved from the LINCS Reverse Search Dashboard\\ref{LINCS Reverse Search Dashboard, https://github.com/MaayanLab/lincs-reverse-search-dashboard}.`,
    legend: `An interactive scatterplot & table showing up and down-regulating perturbation signatures for ${props.inputs?.gene ? props.inputs.gene : 'the gene'} using the RNA-seq-like LINCS L1000 Expression Profile Signatures\\ref{doi:10.1186/s12859-022-04895-5}.`,
  }))
  .build()

export const LINCSL1000ReverseSearchExtractDrugUp = MetaNode('LINCSL1000ReverseSearchExtract[Drug, Up]')
  .meta({
    label: `Extract Up Regulating Perturbagens`,
    description: 'Identify RNA-seq-like LINCS L1000 Chemical Perturbagen Signatures which reverse the expression of the gene.',
    icon: [up_icon],
    external: true,
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
    abstract: `Drugs which up-regulate the expression of ${props.inputs?.search ? props.inputs.search.gene : 'the gene'} were identified from the RNA-seq-like LINCS L1000 Chemical Perturbagens\\ref{doi:10.1093/nar/gkac328}.`,
    introduction: `The LINCS L1000 data\\ref{doi:10.1016/j.cell.2017.10.049} now consists of over 3 million chemical and genetic perturbational signatures. A modified CycleGAN model was used to transorm the L1000 assay measurements to RNA-seq-like expression\\ref{doi:10.1186/s12859-022-04895-5}. Characteristic direction\\ref{doi:10.1186/1471-2105-15-79} was used to compute differential gene expression signatures.`,
    methods: `Pre-computed characteristic direction\\ref{doi:10.1186/1471-2105-15-79}-based differential gene expression signatures from RNA-seq-like LINCS L1000 Expression Profiles\\ref{doi:10.1186/s12859-022-04895-5} which maximally increase or decrease the expression of the gene are retrieved from the LINCS Reverse Search Dashboard\\ref{LINCS Reverse Search Dashboard, https://github.com/MaayanLab/lincs-reverse-search-dashboard}.`,
    legend: `A table of drugs which up-regulate the expression of ${props.inputs?.search ? props.inputs.search.gene : 'the gene'} identified from the RNA-seq-like LINCS L1000 Chemical Perturbagen signatures\\ref{doi:10.1093/nar/gkac328}.`,
  }))
  .build()

export const LINCSL1000ReverseSearchExtractDrugDown = MetaNode('LINCSL1000ReverseSearchExtract[Drug, Down]')
  .meta({
    label: `Extract Down Regulating Perturbagens`,
    description: 'Identify RNA-seq-like LINCS L1000 Chemical Perturbagen Signatures which reverse the expression of the gene.',
    icon: [down_icon],
    external: true,
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
    abstract: `Drugs which down-regulate the expression of ${props.inputs?.search ? props.inputs.search.gene : 'the gene'} were identified from the RNA-seq-like LINCS L1000 Chemical Perturbagens\\ref{doi:10.1093/nar/gkac328}.`,
    introduction: `The LINCS L1000 data\\ref{doi:10.1016/j.cell.2017.10.049} now consists of over 3 million chemical and genetic perturbational signatures. A modified CycleGAN model was used to transorm the L1000 assay measurements to RNA-seq-like expression\\ref{doi:10.1186/s12859-022-04895-5}. Characteristic direction\\ref{doi:10.1186/1471-2105-15-79} was used to compute differential gene expression signatures.`,
    methods: `Pre-computed characteristic direction\\ref{doi:10.1186/1471-2105-15-79}-based differential gene expression signatures from RNA-seq-like LINCS L1000 Expression Profiles\\ref{doi:10.1186/s12859-022-04895-5} which maximally increase or decrease the expression of the gene are retrieved from the LINCS Reverse Search Dashboard\\ref{LINCS Reverse Search Dashboard, https://github.com/MaayanLab/lincs-reverse-search-dashboard}.`,
    legend: `A table of drugs which down-regulate the expression of ${props.inputs?.search ? props.inputs.search.gene : 'the gene'} identified from the RNA-seq-like LINCS L1000 Chemical Perturbagen signatures\\ref{doi:10.1093/nar/gkac328}.`,
  }))
  .build()

export const LINCSL1000ReverseSearchExtractGeneUp = MetaNode('LINCSL1000ReverseSearchExtract[Gene, Up]')
  .meta({
    label: `Extract Up Regulating CRISPR KO Genes`,
    description: 'Identify RNA-seq-like LINCS L1000 CRISPR KO Signatures which reverse the expression of the gene.',
    icon: [up_icon],
    external: true,
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
    abstract: `CRISPR KO genes which up-regulate the expression of ${props.inputs?.search ? props.inputs.search.gene : 'the gene'} were identified from the RNA-seq-like LINCS L1000 CRISPR Knockouts\\ref{doi:10.1093/nar/gkac328}.`,
    introduction: `The LINCS L1000 data\\ref{doi:10.1016/j.cell.2017.10.049} now consists of over 3 million chemical and genetic perturbational signatures. A modified CycleGAN model was used to transorm the L1000 assay measurements to RNA-seq-like expression\\ref{doi:10.1186/s12859-022-04895-5}. Characteristic direction\\ref{doi:10.1186/1471-2105-15-79} was used to compute differential gene expression signatures.`,
    methods: `Pre-computed characteristic direction\\ref{doi:10.1186/1471-2105-15-79}-based differential gene expression signatures from RNA-seq-like LINCS L1000 CRISPR Knockouts\\ref{doi:10.1186/s12859-022-04895-5} which maximally increase or decrease the expression of the gene are retrieved from the LINCS Reverse Search Dashboard\\ref{LINCS Reverse Search Dashboard, https://github.com/MaayanLab/lincs-reverse-search-dashboard}.`,
    legend: `A table of CRISPR KO genes which up-regulate the expression of  ${props.inputs?.search ? props.inputs.search.gene : 'the gene'} identified from the RNA-seq-like LINCS L1000 Chemical Perturbagen signatures\\ref{doi:10.1093/nar/gkac328}.`,
  }))
  .build()

export const LINCSL1000ReverseSearchExtractGeneDown = MetaNode('LINCSL1000ReverseSearchExtract[Gene, Down]')
  .meta({
    label: `Extract Down Regulating CRISPR KOs`,
    description: 'Identify RNA-seq-like LINCS L1000 CRISPR KO Signatures which reverse the expression of the gene.',
    icon: [down_icon],
    external: true,
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
    abstract: `Genes which down-regulate the expression of ${props.inputs?.search ? props.inputs.search.gene : 'the gene'} were identified from the RNA-seq-like LINCS L1000 CRISPR Knockouts\\ref{doi:10.1093/nar/gkac328}.`,
    introduction: `The LINCS L1000 data\\ref{doi:10.1016/j.cell.2017.10.049} now consists of over 3 million chemical and genetic perturbational signatures. A modified CycleGAN model was used to transorm the L1000 assay measurements to RNA-seq-like expression\\ref{doi:10.1186/s12859-022-04895-5}. Characteristic direction\\ref{doi:10.1186/1471-2105-15-79} was used to compute differential gene expression signatures.`,
    methods: `Pre-computed characteristic direction\\ref{doi:10.1186/1471-2105-15-79}-based differential gene expression signatures from RNA-seq-like LINCS L1000 CRISPR Knockouts\\ref{doi:10.1186/s12859-022-04895-5} which maximally increase or decrease the expression of the gene are retrieved from the LINCS Reverse Search Dashboard\\ref{LINCS Reverse Search Dashboard, https://github.com/MaayanLab/lincs-reverse-search-dashboard}.`,
    legend: `A table of CRISPR KO genes which down-regulate the expression of  ${props.inputs?.search ? props.inputs.search.gene : 'the gene'} identified from the RNA-seq-like LINCS L1000 Chemical Perturbagen signatures\\ref{doi:10.1093/nar/gkac328}.`,
  }))
  .build()
