import { GeneCountMatrix, DrugTerm, GeneSet, GeneTerm, PhenotypeTerm, ScoredTissues, Supervenn, ScoredDrugs, ScoredGenes, MetgeneMetabolites, /*VariantSet,*/ LINCSL1000ReverseSearchDashboard, MetGeneSummary, GlyGenResponseNode, MetGeneRxns } from '@/components'

const playbooks = [
  {
    id: '1d66906d-980b-b690-111e-c341145aa080',
    label: 'Use Case 13: Prioritizing Targets for Individual Cancer patients',
    description: 'Given RNA-seq samples from a patient tumor, we screen for targets which are highly expressed in the tumor but lowly expressed across all healthy human tissues in GTEx. Detailed information about the the selected target are queried from several DCCs.',
    gpt_summary: 'The process began with the upload of a file, which was then parsed as an input gene count matrix. The next step involved identifying significantly over-expressed genes in comparison to normal tissue in GTEx, resulting in the selection of IMP3 for further investigation.\n\nTo identify drugs that down-regulate IMP3 expression, RNA-seq-like LINCS L1000 Chemical Perturbagens were utilized. Additionally, genes that down-regulate IMP3 were identified from RNA-seq-like LINCS L1000 CRISPR Knockouts. The list of genes was then filtered by IDG Understudied Proteins.\n\nThe Metabolomics Workbench was used to identify associated metabolites and relevant reactions for IMP3. Regulatory elements were obtained from the Linked Data Hub. Lastly, the GlyGen database was searched to identify a relevant set of proteins originating from IMP3.',
    published: 'Mar 30, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CB-BY-NC-SA-4.0',
    url: 'https://playbook-workflow-builder.cloud/report/1d66906d-980b-b690-111e-c341145aa080',
    dataSources: [
      'KidsFirst',
      'LINCS L1000',
      'GTEx',
      'exRNA',
      'GlyGen',
      'Metabolomics',
      'IDG',
    ],
    inputs: [
      GeneCountMatrix,
    ],
    outputs: [
      LINCSL1000ReverseSearchDashboard,
      MetGeneSummary,
      ScoredGenes,
      ScoredDrugs,
      MetGeneRxns,
      MetgeneMetabolites,
      /*VariantSet,*/
      GlyGenResponseNode,
    ],
    collapsed: {
      4: false,
      5: false,
      7: false,
      8: false,
      9: false,
      10: false,
      11: false,
      12: false,
    } as Record<number, boolean>,
    clicks: 5,
  },
  {
    id: '0db222e1-9958-b01d-4e22-fef5599ce1f5',
    label: 'Use Case 4: Identify the Tissue Activity for a TF based on its Targets',
    description: 'Given a Transcription Factor, we collect its targets by various resources and then Enrich the set of consensus targets against GTEx Tissue expression. The result is ranked tissues potentially regulated by the transcription factor.',
    gpt_summary: 'To start the workflow, KLF4 is selected as the search term. Matching terms from three libraries are identified: ENCODE TF ChIP-seq 2015, ChEA 2022, and ARCHS4 TF Co-Expression. The matching terms are assembled into gene sets, which are then combined into one library. A consensus gene set is created by retaining genes that appear in at least two sets. The gene set is enriched against the GTEx Tissues V8 2023 library to identify statistically significant GTEx Tissue Signatures. The libraries used in this workflow are: \n1. An integrated encyclopedia of DNA elements in the human genome. Nature vol. 489 57–74 (2012). doi:10.1038/nature11247 \n2. Keenan, A. B. et al. ChEA3: transcription factor enrichment analysis by orthogonal omics integration. Nucleic Acids Research vol. 47 W212–W224 (2019). doi:10.1093/nar/gkz446 \n3. Lachmann, A. et al. Massive mining of publicly available RNA-seq data from human and mouse. Nature Communications vol. 9 (2018). doi:10.1038/s41467-018-03751-6 \n4. Lonsdale, J. et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics vol. 45 580–585 (2013). doi:10.1038/ng.2653',
    published: 'Mar 16, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CB-BY-NC-SA-4.0',
    url: 'https://playbook-workflow-builder.cloud/report/0db222e1-9958-b01d-4e22-fef5599ce1f5',
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
    collapsed: {} as Record<number, boolean>,
    clicks: 2,
  },
  {
    id: '0c0e357d-79b7-47d0-ec94-7915a88bf493',
    label: 'Use Case 2: Explain MOAs of Side Effects for Approved Drugs',
    description: 'For a side effect and a drug, I would like to know if there are genes from the LINCS L1000 experiments that are up or down-regulated by the drugs that are also known to be involved with the side effect based on literature co-mentions or GWAS. I would like to know if such overlap is statistically significant. I would also like the results to be visualized using a SuperVenn diagram.',
    gpt_summary: 'To start the workflow, the search term \"atrial fibrillation\" was selected. The matching terms from the MGI Mammalian Phenotype Level 4 2019 library and GWAS Catalog 2019 library were identified and assembled into gene sets. A consensus gene set was created by only keeping genes that appeared in at least two sets.\n\nSimilarly, for the search term \"Ibrutinib,\" the matching terms from the LINCS L1000 Chem Pert Consensus Sigs library were identified and assembled into gene sets. All the identified gene sets were combined into one gene set library.\n\nFinally, a Supervenn diagram was used to visualize the collection of gene sets. The sources for the libraries used in this workflow are: \n1. Mouse Genome Database (MGD)\n2. NHGRI-EBI GWAS Catalog\n3. SigCom LINCS.',
    published: 'Mar 16, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CB-BY-NC-SA-4.0',
    url: 'https://playbook-workflow-builder.cloud/report/0c0e357d-79b7-47d0-ec94-7915a88bf493',
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
    collapsed: {} as Record<number, boolean>,
    clicks: 1,
  },
]

export default playbooks
