import { GeneCountMatrix, DrugTerm, GeneSet, GeneTerm, PhenotypeTerm, ScoredTissues, Supervenn, ScoredDrugs, ScoredGenes, ProteinProductInformation, MetgeneMetabolites, VariantSet } from '@/components'

const playbooks = [
  {
    id: '9c14c4b4-a7f4-4433-ee30-b9a3abb6439a',
    label: 'Use Case 13: Prioritizing Targets for Individual Cancer patients',
    description: 'Given RNA-seq samples from a patient tumor, we screen for targets which are highly expressed in the tumor but lowly expressed across all healthy human tissues in GTEx. Detailed information about the the selected target are queried from several DCCs.',
    gpt_summary: 'A file containing samples from a lung cancer patient was uploaded and parsed as an input gene count matrix. Significantly over-expressed genes when compared to normal tissue in GTEx were identified, and IMP3 was chosen for further investigation. Genes that up-regulate the expression of IMP3 were identified from RNA-seq-like LINCS L1000 CRISPR Knockouts [1], while drugs that down-regulate the expression of IMP3 were identified from RNA-seq-like LINCS L1000 Chemical Perturbagens [1]. More information about the gene was obtained with the MyGene.info API [2,3], and regulatory elements were obtained from the Linked Data Hub [4]. Next, the GlyGen database [5] was searched to identify a relevant set of proteins that originate from the gene. Finally, the gene was searched in the Metabolomics Workbench [6] to identify associated metabolites.\n\nReferences:\n1. Evangelista, J. E. et al. SigCom LINCS: data and metadata search engine for a million gene expression signatures. Nucleic Acids Research vol. 50 W697–W709 (2022). doi:10.1093/nar/gkac328\n2. Xin, J. et al. High-performance web services for querying gene and variant annotation. Genome Biology vol. 17 (2016). doi:10.1186/s13059-016-0953-9\n3. Wu, C., MacLeod, I. & Su, A. I. BioGPS and MyGene.info: organizing online, gene-centric information. Nucleic Acids Research vol. 41 D561–D565 (2012). doi:10.1093/nar/gks1114\n4. Linked Data Hub, https://genboree.org/cfde-gene-dev/\n5. York, W. S. et al. GlyGen: Computational and Informatics Resources for Glycoscience. Glycobiology vol. 30 72–73 (2019). doi:10.1093/glycob/cwz080\n6. The Metabolomics Workbench, https://www.metabolomicsworkbench.org/',
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
    gpt_summary: 'To start the workflow, the search term KLF4 is selected. Matching terms from three different libraries are identified: the ENCODE TF ChIP-seq 2015 library, the ChEA 2022 library, and the ARCHS4 TF Co-Expression library. The identified matching terms from each library are collected into separate gene sets. These gene sets are combined to create one gene set library. A consensus gene set is then formed by only retaining genes that appear in at least two sets. Finally, the gene set is enriched against the GTEx Tissues V8 2023 library to identify statistically significant GTEx Tissue Signatures. The sources for the libraries used in this workflow are: \n\n1. http://genome.ucsc.edu/ENCODE/downloads.html\n2. Keenan, A. B. et al. ChEA3: transcription factor enrichment analysis by orthogonal omics integration. Nucleic Acids Research vol. 47 W212–W224 (2019). doi:10.1093/nar/gkz446\n3. Lachmann, A. et al. Massive mining of publicly available RNA-seq data from human and mouse. Nature Communications vol. 9 (2018). doi:10.1038/s41467-018-03751-6\n4. Lonsdale, J. et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics vol. 45 580–585 (2013). doi:10.1038/ng.2653',
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
    gpt_summary: 'To start the workflow, the search term \"atrial fibrillation\" is selected. Matching terms from the MGI Mammalian Phenotype Level 4 2019 library and the GWAS Catalog 2019 library are identified and assembled into gene sets. A consensus gene set is created by retaining only genes that appear in at least two sets.\n\nSimilarly, for the search term \"Ibrutinib\", matching terms from the LINCS L1000 Chem Pert Consensus Sigs library are assembled into gene sets. The collected gene sets are then combined into one gene set library.\n\nFinally, the collection of gene sets is visualized using a Supervenn diagram. All the necessary sources for this workflow are cited in the text.',
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
