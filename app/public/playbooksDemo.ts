import { GeneCountMatrix, DrugTerm, GeneSet, GeneTerm, PhenotypeTerm, ScoredTissues, Supervenn, ScoredDrugs, ScoredGenes, ProteinProductInformation, MetgeneMetabolites, VariantSet } from '@/components'

const playbooks = [
  {
    id: '9c14c4b4-a7f4-4433-ee30-b9a3abb6439a',
    label: 'Use Case 13: Prioritizing Targets for Individual Cancer patients',
    description: 'Given RNA-seq samples from a patient tumor, we screen for targets which are highly expressed in the tumor but lowly expressed across all healthy human tissues in GTEx. Detailed information about the the selected target are queried from several DCCs.',
    gpt_summary: undefined,
    published: 'Mar 27, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CB-BY-NC-SA-4.0',
    url: 'https://dev.playbook-workflow-builder.cloud/report/9c14c4b4-a7f4-4433-ee30-b9a3abb6439a',
    dataSources: [
      'KidsFirst',
      'LINCS L1000',
      'GTEx',
      'exRNA',
      'GlyGen',
      'Metabolomics',
    ],
    inputs: [
      GeneCountMatrix,
    ],
    outputs: [
      ScoredDrugs,
      ScoredGenes,
      VariantSet,
      MetgeneMetabolites,
      ProteinProductInformation,
    ],
    clicks: 5,
  },
  {
    id: '0db222e1-9958-b01d-4e22-fef5599ce1f5',
    label: 'Use Case 4: Identify the Tissue Activity for a TF based on its Targets',
    description: 'Given a Transcription Factor, we collect its targets by various resources and then Enrich the set of consensus targets against GTEx Tissue expression. The result is ranked tissues potentially regulated by the transcription factor.',
    gpt_summary: undefined,
    published: 'Mar 16, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CB-BY-NC-SA-4.0',
    url: 'https://dev.playbook-workflow-builder.cloud/report/0db222e1-9958-b01d-4e22-fef5599ce1f5',
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
    id: '0c0e357d-79b7-47d0-ec94-7915a88bf493',
    label: 'Use Case 2: Explain MOAs of Side Effects for Approved Drugs',
    description: 'For a side effect and a drug, I would like to know if there are genes from the LINCS L1000 experiments that are up or down-regulated by the drugs that are also known to be involved with the side effect based on literature co-mentions or GWAS. I would like to know if such overlap is statistically significant. I would also like the results to be visualized using a SuperVenn diagram.',
    gpt_summary: undefined,
    published: 'Mar 16, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CB-BY-NC-SA-4.0',
    url: 'https://dev.playbook-workflow-builder.cloud/report/0c0e357d-79b7-47d0-ec94-7915a88bf493',
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
