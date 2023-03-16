import { DrugTerm, GeneSet, GeneTerm, PhenotypeTerm, ScoredTissues, Supervenn } from '@/components'

const playbooks = [
  {
    id: '9abb4a04-8bcd-4ecf-1628-6c9821e16ccc',
    label: 'Use Case 4: Identify the Tissue Activity for a TF based on its Targets',
    description: 'Given a Transcription Factor, we collect its targets by various resources and then Enrich the set of consensus targets against GTEx Tissue expression. The result is ranked tissues potentially regulated by the transcription factor.',
    published: 'Mar 16, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CB-BY-NC-SA-4.0',
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
    clicks: 2,
  },
  {
    id: 'd30f6646-f29b-4427-41bb-641be48805a4',
    label: 'Use Case 2: Explain MOAs of Side Effects for Approved Drugs',
    description: 'For a side effect and a drug, I would like to know if there are genes from the LINCS L1000 experiments that are up or down-regulated by the drugs that are also known to be involved with the side effect based on literature co-mentions or GWAS. I would like to know if such overlap is statistically significant. I would also like the results to be visualized using a SuperVenn diagram.',
    published: 'Mar 16, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CB-BY-NC-SA-4.0',
    url: 'https://playbook-partnership.dev.maayanlab.cloud/graph/c28a9d39-077f-43d7-ad9b-af961a3055d7',
    dataSources: [
      'GWAS Catalog',
      'KOMP',
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
    clicks: 1,
  },
]

export default playbooks
