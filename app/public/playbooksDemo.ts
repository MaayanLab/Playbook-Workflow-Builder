import { DrugTerm, GeneSet, GeneTerm, PhenotypeTerm, ScoredTissues, Supervenn } from '@/components'

const playbooks = [
  {
    id: '9abb4a04-8bcd-4ecf-1628-6c9821e16ccc',
    label: 'GTEx Tissue Enrichment of Transcription Factor Targets',
    description: 'Given a Transcription Factor, we collect its targets by various resources and then Enrich the set of consensus targets against GTEx Tissue expression. The result is ranked tissues potentially regulated by the transcription factor.',
    url: 'https://playbook-partnership.maayanlab.cloud/report/9abb4a04-8bcd-4ecf-1628-6c9821e16ccc',
    dataSources: [
      'ENCODE',
      'ChEA',
      'ARCHS4',
      'GTEx',
    ],
    inputs: [
      GeneTerm,
    ],
    outputs: [
      ScoredTissues,
    ],
  },
  {
    id: 'd30f6646-f29b-4427-41bb-641be48805a4',
    label: 'Explain MOAs of Side Effects for Approved Drugs',
    description: 'For a side effect and a drug, I would like to know if there are genes from the LINCS L1000 experiments that are up or down-regulated by the drugs that are also known to be involved with the side effect based on literature co-mentions or GWAS. I would like to know if such overlap is statistically significant. I would also like the results to be visualized using a SuperVenn diagram.',
    url: 'https://playbook-partnership.maayanlab.cloud/report/d30f6646-f29b-4427-41bb-641be48805a4',
    dataSources: [
      'GWAS Catalog',
      'GeneRif',
      'LINCS L1000',
    ],
    inputs: [
      PhenotypeTerm,
      DrugTerm,
    ],
    outputs: [
      GeneSet,
      Supervenn,
    ],
  },
]
export default playbooks