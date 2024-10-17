import { MetgeneMetaboliteTable } from '@/components/MW/metgene_metabolite_table'
import { MetGeneRxnTable } from '@/components/MW/metgene_rxn_table'
import { MetGeneStudyTable } from '@/components/MW/metgene_study_table'
import { MetGeneSummary } from '@/components/MW/metgene_summary'
import { ScoredDrugs, ScoredGenes, ScoredTissues } from '@/components/core/scored'
import { GeneSet, RegulatoryElementSet, VariantSet } from '@/components/core/set'
import { DrugTerm, GeneTerm, PhenotypeTerm, VariantTerm } from '@/components/core/term'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { GeneSignature } from '@/components/data/gene_signature'
import { GlyGenProteinResponseNode } from '@/components/gly_gen'
import { LINCSL1000ReverseSearchDashboard } from '@/components/lincs/l1000-reverse-search'
import { Supervenn } from '@/components/service/hyposet'
import { PlotlyPlot } from '@/components/viz/plotly'

const playbooks = [
  {
    id: '6a7af61b-ef0a-7687-5d6e-02deeb253172',
    label: 'Use Case 1: Explain Drug-Drug Interactions',
    description: `Give two drugs and an adverse event that is known to be caused by the drug-drug interactions, I would like to know if there are overlapping genes between genes that are either up or down regulated by the drugs from LINCS and genes associated with the adverse event either based on GWAS, gene mentions in the literature, and genes associated with mouse and human phenotypes. I would like to also know if the overlap between these genes is statistically significant. I would also like to have the results visualized as a Venn diagram.`,
    gpt_summary: `To start the workflow, the search term "Inflammation" is selected. This is followed by selecting the search term "Penicillin" and then "Cortisol".

From the Enrichr database, gene sets with set labels containing "Inflammation" are queried. The identified matching terms from the GWAS Catalog 2019 library are assembled into a collection of gene sets. A GMT (gene matrix transposed) is extracted from the Enrichr results for GWAS_Catalog_2019. Then, all the identified gene sets are combined using the union set operation.

Similarly, gene sets with set labels containing "Penicillin" are queried from Enrichr. The identified matching terms from the LINCS L1000 Chem Pert Consensus Sigs library are assembled into a collection of gene sets. A GMT is extracted from the Enrichr results for LINCS_L1000_Chem_Pert_Consensus_Sigs.

Next, gene sets with set labels containing "Cortisol" are queried from Enrichr. The identified matching terms from the LINCS L1000 Chem Pert Consensus Sigs library are assembled into a collection of gene sets. A GMT is extracted from the Enrichr results for LINCS_L1000_Chem_Pert_Consensus_Sigs.

All the identified gene sets from the three search terms are combined into one gene set library. This collection of gene sets is then visualized with a Supervenn diagram.

Here are the references for the sources used:
1. Xie, Z. et al. Gene Set Knowledge Discovery with Enrichr. Current Protocols vol. 1 (2021). doi:10.1002/cpz1.90
2. Sollis, E. et al. The NHGRI-EBI GWAS Catalog: knowledgebase and deposition resource. Nucleic Acids Research vol. 51 D977–D985 (2022). doi:10.1093/nar/gkac1010
3. Blake, J. A. et al. Mouse Genome Database (MGD): Knowledgebase for mouse–human comparative biology. Nucleic Acids Research vol. 49 D981–D987 (2020). doi:10.1093/nar/gkaa1083
4. Köhler, S. et al. The Human Phenotype Ontology in 2021. Nucleic Acids Research vol. 49 D1207–D1217 (2020). doi:10.1093/nar/gkaa1043
5. Evangelista, J. E. et al. SigCom LINCS: data and metadata search engine for a million gene expression signatures. Nucleic Acids Research vol. 50 W697–W709 (2022). doi:10.1093/nar/gkac328`,
    published: 'July 26, 2023',
    version: '1.0.1',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CC-BY-NC-SA-4.0',
    dataSources: [
      'KOMP',
      'HPO',
      'GWAS Catalog',
      'Geneshot',
      'LINCS',
    ],
    inputs: [
      PhenotypeTerm,
      DrugTerm,
      DrugTerm,
    ],
    outputs: [
      Supervenn,
    ],
    workflow: {"data":{"e34973a6-553e-caa9-2f19-44efaf82557e":{"type":"Term[Phenotype]","value":"Inflammation"},"909e6678-c8ed-ec49-4577-2911d9ebedff":{"type":"Term[Drug]","value":"Penicillin"},"77232520-9ad0-83a6-49c6-ca5c3f774ad5":{"type":"Term[Drug]","value":"Cortisol"},"ce129a35-fb18-b067-ba63-21baf5eb44df":{"type":"GenesetsToGMT","value":{"terms":{"0":"GWAS","1":"MGI","2":"HPO"},"descriptions":{"0":"","1":"","2":""}}}},"workflow":[{"id":"48fcbe58-ca9c-a9bc-2785-76ffeb1f4f7d","type":"Input[Phenotype]","data":{"id":"e34973a6-553e-caa9-2f19-44efaf82557e"}},{"id":"b79e3919-1d40-a835-bbe0-b5c226b9e2ff","type":"Input[Drug]","data":{"id":"909e6678-c8ed-ec49-4577-2911d9ebedff"}},{"id":"073dd744-9101-d2cb-28f5-c08e60349f95","type":"Input[Drug]","data":{"id":"77232520-9ad0-83a6-49c6-ca5c3f774ad5"}},{"id":"4dd140a2-9fd3-de1b-4b3b-ea79315335a3","type":"EnrichrTermSearch[Phenotype]","inputs":{"term":{"id":"48fcbe58-ca9c-a9bc-2785-76ffeb1f4f7d"}}},{"id":"ceefad3e-6395-921e-3870-fecd1d077b0c","type":"ExtractEnrichrTermSearch[GWAS_Catalog_2019]","inputs":{"searchResults":{"id":"4dd140a2-9fd3-de1b-4b3b-ea79315335a3"}}},{"id":"316aed3b-e820-a959-93b6-44c8dd2b14fa","type":"EnrichrSetTToGMT[Phenotype]","inputs":{"enrichrset":{"id":"ceefad3e-6395-921e-3870-fecd1d077b0c"}}},{"id":"c92d494e-b9c2-4cb8-eeaa-30b8c857edff","type":"GMTUnion","inputs":{"gmt":{"id":"316aed3b-e820-a959-93b6-44c8dd2b14fa"}}},{"id":"ad12f973-227c-b5cc-b185-7dc9fad7f102","type":"ExtractEnrichrTermSearch[MGI_Mammalian_Phenotype_Level_4_2019]","inputs":{"searchResults":{"id":"4dd140a2-9fd3-de1b-4b3b-ea79315335a3"}}},{"id":"dd0e5963-b9f8-59d5-df4b-e3b08d1c00e7","type":"EnrichrSetTToGMT[Phenotype]","inputs":{"enrichrset":{"id":"ad12f973-227c-b5cc-b185-7dc9fad7f102"}}},{"id":"be2277af-b1d5-4c3f-881e-3de4bbade0e6","type":"GMTUnion","inputs":{"gmt":{"id":"dd0e5963-b9f8-59d5-df4b-e3b08d1c00e7"}}},{"id":"1a3378d9-e5f5-3c81-eb37-a7c0d2f83192","type":"ExtractEnrichrTermSearch[Human_Phenotype_Ontology]","inputs":{"searchResults":{"id":"4dd140a2-9fd3-de1b-4b3b-ea79315335a3"}}},{"id":"b705e386-e46f-a6d2-864f-8cfa1ca614fe","type":"EnrichrSetTToGMT[Phenotype]","inputs":{"enrichrset":{"id":"1a3378d9-e5f5-3c81-eb37-a7c0d2f83192"}}},{"id":"674a5cea-5334-5767-49ff-60f1f04ab54a","type":"GMTUnion","inputs":{"gmt":{"id":"b705e386-e46f-a6d2-864f-8cfa1ca614fe"}}},{"id":"23365a0b-4869-9c36-3a9c-abc9b3fd97c0","type":"EnrichrTermSearch[Drug]","inputs":{"term":{"id":"b79e3919-1d40-a835-bbe0-b5c226b9e2ff"}}},{"id":"113875b7-f3b2-413e-68e0-c7732dae6523","type":"ExtractEnrichrTermSearch[LINCS_L1000_Chem_Pert_Consensus_Sigs]","inputs":{"searchResults":{"id":"23365a0b-4869-9c36-3a9c-abc9b3fd97c0"}}},{"id":"7b60355d-9b2e-034b-517c-d182a595a88b","type":"EnrichrSetTToGMT[Drug]","inputs":{"enrichrset":{"id":"113875b7-f3b2-413e-68e0-c7732dae6523"}}},{"id":"d1ba51b9-568d-954b-25ef-c871cebc2e5c","type":"EnrichrTermSearch[Drug]","inputs":{"term":{"id":"073dd744-9101-d2cb-28f5-c08e60349f95"}}},{"id":"68878725-8b82-90eb-f93d-6d4814470ec9","type":"ExtractEnrichrTermSearch[LINCS_L1000_Chem_Pert_Consensus_Sigs]","inputs":{"searchResults":{"id":"d1ba51b9-568d-954b-25ef-c871cebc2e5c"}}},{"id":"e6ccb2b4-d7f4-ecb3-a536-28fdc2845b85","type":"EnrichrSetTToGMT[Drug]","inputs":{"enrichrset":{"id":"68878725-8b82-90eb-f93d-6d4814470ec9"}}},{"id":"edb44eff-e405-8827-2d6b-fa4c021cb352","type":"GenesetsToGMT","inputs":{"genesets:0":{"id":"c92d494e-b9c2-4cb8-eeaa-30b8c857edff"},"genesets:1":{"id":"be2277af-b1d5-4c3f-881e-3de4bbade0e6"},"genesets:2":{"id":"674a5cea-5334-5767-49ff-60f1f04ab54a"}},"data":{"id":"ce129a35-fb18-b067-ba63-21baf5eb44df"}},{"id":"f4df6163-5723-6acc-93dd-4f3cc719d3d1","type":"GMTConcatenate","inputs":{"gmts:0":{"id":"7b60355d-9b2e-034b-517c-d182a595a88b"},"gmts:1":{"id":"e6ccb2b4-d7f4-ecb3-a536-28fdc2845b85"},"gmts:2":{"id":"edb44eff-e405-8827-2d6b-fa4c021cb352"}}},{"id":"80102c29-34f5-2fb2-bf7f-650ac365d518","type":"SupervennFromGMT","inputs":{"gmt":{"id":"f4df6163-5723-6acc-93dd-4f3cc719d3d1"}}}],"metadata":{"id":"fcce635a-fd1f-76f8-7fd2-5d909b43fbf4","title":"Use Case 1: Explain Drug-Drug Interactions","description":"Give two drugs and an adverse event that is known to be caused by the drug-drug interactions, I would like to know if there are overlapping genes between genes that are either up or down regulated by the drugs from LINCS and genes associated with the adverse event either based on GWAS, gene mentions in the literature, and genes associated with mouse and human phenotypes. I would like to also know if the overlap between these genes is statistically significant. I would also like to have the results visualized as a Venn diagram.","summary":"auto","gpt_summary":""}},
  },
  {
    id: 'c4d40504-57b6-d48f-6d12-d47891e26f2d',
    label: 'Use Case 2: Explain MOAs of Side Effects for Approved Drugs',
    description: `For a side effect and a drug, find differentially expressed genes from the LINCS L1000 resource that are up- or down-regulated by the drug, and are also associated with the side-effect based on literature co-mentions or GWAS. If overlapping genes are found, compute whether such overlap is statistically significant and visualize the results with a supervenn diagram.`,
    gpt_summary: `To start the workflow, the search term \"atrial fibrillation\" was selected. The matching terms from the MGI Mammalian Phenotype Level 4 2019 library and GWAS Catalog 2019 library were identified and assembled into gene sets. A consensus gene set was created by only keeping genes that appeared in at least two sets.

Similarly, for the search term \"Ibrutinib,\" the matching terms from the LINCS L1000 Chem Pert Consensus Sigs library were identified and assembled into gene sets. All the identified gene sets were combined into one gene set library.

Finally, a Supervenn diagram was used to visualize the collection of gene sets. The sources for the libraries used in this workflow are: 
1. Mouse Genome Database (MGD)
2. NHGRI-EBI GWAS Catalog
3. SigCom LINCS.`,
    published: 'Mar 16, 2023',
    version: '1.0.1',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CC-BY-NC-SA-4.0',
    dataSources: [
      'GWAS Catalog',
      'KOMP',
      'LINCS',
    ],
    inputs: [
      PhenotypeTerm,
      DrugTerm,
    ],
    outputs: [
      GeneSet,
      Supervenn,
    ],
    workflow: {"data":{"3789e8b8-6190-f7d3-5289-2200cea1ce26":{"type":"Input[Disease]","value":"atrial fibrillation"},"cc64cfd5-d920-0e0a-45e4-5e848efd0db9":{"type":"Input[Drug]","value":"Ibrutinib"},"bf6292f6-ea5e-8d04-3330-e289d2826a41":{"type":"GenesetsToGMT","value":{"terms":{"0":"MGI consensus AF","1":"GWAS Catalog consensus AF"},"descriptions":{"0":"","1":""}}}},"workflow":[{"id":"80d31c4f-c14c-d091-ba6f-535b0e67a09f","type":"Input[Disease]","data":{"id":"3789e8b8-6190-f7d3-5289-2200cea1ce26"}},{"id":"1f2918be-67e5-ff5d-9d31-05babd9fcebe","type":"Input[Drug]","data":{"id":"cc64cfd5-d920-0e0a-45e4-5e848efd0db9"}},{"id":"d2df26ee-3804-c797-cb54-6ab8fa25a452","type":"EnrichrTermSearch[Disease]","inputs":{"term":{"id":"80d31c4f-c14c-d091-ba6f-535b0e67a09f"}}},{"id":"03682246-d326-ba49-e440-30dd6a44e6f0","type":"ExtractEnrichrTermSearch[MGI_Mammalian_Phenotype_Level_4_2021]","inputs":{"searchResults":{"id":"d2df26ee-3804-c797-cb54-6ab8fa25a452"}}},{"id":"fa24037f-20ab-6d01-f7f4-abbdb5a93374","type":"EnrichrSetTToGMT[Phenotype]","inputs":{"enrichrset":{"id":"03682246-d326-ba49-e440-30dd6a44e6f0"}}},{"id":"53f2f91c-78c9-509b-ecc6-e492b7bd1692","type":"GMTConsensus","inputs":{"gmt":{"id":"fa24037f-20ab-6d01-f7f4-abbdb5a93374"}}},{"id":"9c061636-7d66-db0d-78bf-9ca7bf5eaed5","type":"ExtractEnrichrTermSearch[GWAS_Catalog_2019]","inputs":{"searchResults":{"id":"d2df26ee-3804-c797-cb54-6ab8fa25a452"}}},{"id":"4ded35b6-88d7-70d8-8f0f-85bd2ca12586","type":"EnrichrSetTToGMT[Phenotype]","inputs":{"enrichrset":{"id":"9c061636-7d66-db0d-78bf-9ca7bf5eaed5"}}},{"id":"bd505905-922c-7cc8-8120-113fb05c9161","type":"GMTConsensus","inputs":{"gmt":{"id":"4ded35b6-88d7-70d8-8f0f-85bd2ca12586"}}},{"id":"c55e1630-7022-49ad-47a7-f40e994034a4","type":"GenesetsToGMT","inputs":{"genesets:0":{"id":"53f2f91c-78c9-509b-ecc6-e492b7bd1692"},"genesets:1":{"id":"bd505905-922c-7cc8-8120-113fb05c9161"}},"data":{"id":"bf6292f6-ea5e-8d04-3330-e289d2826a41"}},{"id":"cb0f4acb-672d-12c8-24a5-3b48942a94aa","type":"EnrichrTermSearch[Drug]","inputs":{"term":{"id":"1f2918be-67e5-ff5d-9d31-05babd9fcebe"}}},{"id":"a632aeb4-bacf-d014-3af1-6c713a970389","type":"ExtractEnrichrTermSearch[LINCS_L1000_Chem_Pert_Consensus_Sigs]","inputs":{"searchResults":{"id":"cb0f4acb-672d-12c8-24a5-3b48942a94aa"}}},{"id":"1fa4e6da-d5be-8054-8683-87bede64a65a","type":"EnrichrSetTToGMT[Drug]","inputs":{"enrichrset":{"id":"a632aeb4-bacf-d014-3af1-6c713a970389"}}},{"id":"6d0e6a1a-956d-5322-85cd-2987ed9b15c4","type":"GMTConcatenate","inputs":{"gmts:0":{"id":"c55e1630-7022-49ad-47a7-f40e994034a4"},"gmts:1":{"id":"1fa4e6da-d5be-8054-8683-87bede64a65a"}}},{"id":"5f4751a5-759f-a9f6-c615-c73662df5170","type":"SupervennFromGMT","inputs":{"gmt":{"id":"6d0e6a1a-956d-5322-85cd-2987ed9b15c4"}}}],"metadata":{"id":"31b447be-b6fc-1e4b-74e9-9208174f2984","title":"Use Case 2: Explain MOAs of Side Effects for Approved Drugs","description":"For a side effect and a drug, find differentially expressed genes from the LINCS L1000 resource that are up- or down-regulated by the drug, and are also associated with the side-effect based on literature co-mentions or GWAS. If overlapping genes are found, compute whether such overlap is statistically significant and visualize the results with a supervenn diagram.","summary":"auto","gpt_summary":""}},
  },
  {
    id: '4153d6c3-0d1f-ba10-eeb6-c7792318c042',
    label: 'Use Case 3: Compounds to Reverse Disease Signatures',
    description: `Using differential expression signatures from two independent sources, identify and rank Consensus L1000 Small Molecules capable of reversing the expression signatures. We consider a case study involving Aging signatures between GTEx and GEO, drugs proposed should be considered as potential treatments for aging.`,
    gpt_summary: `The following steps were taken:

1. The GEO Aging Signatures file was uploaded and loaded as a gene signature.
2. The GTEx Aging Signatures file was uploaded and loaded as a gene signature.
3. Significant genes were extracted from the GEO Aging Signatures.
4. Significant genes were extracted from the GTEx Aging Signatures.
5. Reversers and mimickers were identified from over 1 million signatures using SigCom LINCS [1].
6. Resolved drugs were obtained from the LINCS L1000 Chemical Perturbagens library.
7. Reversers and mimickers were identified from over 1 million signatures using SigCom LINCS [1].
8. Mean scores were computed.
9. The drugs were filtered by FDA Approved Drugs using the PubChem APIs [2].

References:
1. Evangelista, J. E. et al. SigCom LINCS: data and metadata search engine for a million gene expression signatures. Nucleic Acids Research vol. 50 W697–W709 (2022). doi:10.1093/nar/gkac328
2. Kim, S. et al. PubChem 2023 update. Nucleic Acids Research vol. 51 D1373–D1380 (2022). doi:10.1093/nar/gkac956`,
    published: 'July 26, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CC-BY-NC-SA-4.0',
    dataSources: [
      'GTEx',
      'GEO',
      'LINCS',
    ],
    inputs: [
      GeneSignature,
      GeneSignature,
    ],
    outputs: [
      ScoredDrugs,
    ],
    workflow: {"data":{"03e58af8-72c8-2ead-3edd-54c5fd6fefd8":{"type":"FileURL","value":{"description":"GEO Aging Signatures","url":"file:///tmp/2mtsGk5TBg8SA4YlKZ00WwU6.tsv","filename":"GEO_aging_signature_limma.tsv","sha256":"1eccf6972a36d79593dd0f09a297a15cf161295372e2501216d4a01ca9c348e2","size":44677}},"c699f27c-6be9-bca8-a11b-97cf89fe4530":{"type":"FileURL","value":{"description":"GTEx Aging Signatures","url":"file:///tmp/WIefGM0HpMHA6jdsFs-pN1ks.tsv","filename":"GTEx_aging_signature_limma.tsv","sha256":"33e522bf1b556ef6cbb2603dd0917003c94fa2b1424a425e5ccf95c397455365","size":43161}}},"workflow":[{"id":"826c33f7-a796-d882-55af-bbfb7180f7fb","type":"FileInput","data":{"id":"03e58af8-72c8-2ead-3edd-54c5fd6fefd8"}},{"id":"71636ce3-191a-ee9a-87e2-a9e15595f9fd","type":"GeneSigFromFile","inputs":{"file":{"id":"826c33f7-a796-d882-55af-bbfb7180f7fb"}}},{"id":"7c3b27a0-231d-ba81-5385-7270aa2030ca","type":"FileInput","data":{"id":"c699f27c-6be9-bca8-a11b-97cf89fe4530"}},{"id":"3ac4c848-0f96-c0ee-d2ac-935b29e93a18","type":"GeneSigFromFile","inputs":{"file":{"id":"7c3b27a0-231d-ba81-5385-7270aa2030ca"}}},{"id":"9c6bcd82-e98b-527e-40a2-26c7cb3ed3f7","type":"ScoredGenesFromSignature","inputs":{"sig":{"id":"71636ce3-191a-ee9a-87e2-a9e15595f9fd"}}},{"id":"08939252-4437-42ef-db51-0abc4898897a","type":"ScoredGenesFromSignature","inputs":{"sig":{"id":"3ac4c848-0f96-c0ee-d2ac-935b29e93a18"}}},{"id":"cff25269-e005-63e8-0cae-d7116b763446","type":"SigComLINCSSignatureSearch","inputs":{"genes":{"id":"08939252-4437-42ef-db51-0abc4898897a"}}},{"id":"108ee6f3-5eb6-5cb8-b626-910f2a4e57a8","type":"ExtractSigComLINCSSignatureSearch[LINCS L1000 Chemical Perturbagens]","inputs":{"searchResults":{"id":"cff25269-e005-63e8-0cae-d7116b763446"}}},{"id":"7ffae81b-116f-2248-5cd8-6cc7485671d5","type":"SigComLINCSSignatureSearch","inputs":{"genes":{"id":"9c6bcd82-e98b-527e-40a2-26c7cb3ed3f7"}}},{"id":"d65edecb-c42b-7d25-4805-f71994875ab7","type":"ExtractSigComLINCSSignatureSearch[LINCS L1000 Chemical Perturbagens]","inputs":{"searchResults":{"id":"7ffae81b-116f-2248-5cd8-6cc7485671d5"}}},{"id":"bcfc5829-b23a-e907-b7db-1d1f02050f5f","type":"MeanScoredT[Drug]","inputs":{"scored:0":{"id":"108ee6f3-5eb6-5cb8-b626-910f2a4e57a8"},"scored:1":{"id":"d65edecb-c42b-7d25-4805-f71994875ab7"}}},{"id":"c341af16-54c2-a92b-fa2c-bf485f368073","type":"FilterFDAApprovedDrugs[Scored[Drug]]","inputs":{"drugs":{"id":"bcfc5829-b23a-e907-b7db-1d1f02050f5f"}}}],"metadata":{"id":"ad8bc25c-afaa-ba7b-b7b2-80035b5da098","title":"Use Case 3: Compounds to Reverse Disease Signatures","description":"Using differential expression signatures from two independent sources, identify and rank Consensus L1000 Small Molecules capable of reversing the expression signatures. We consider a case study involving Aging signatures between GTEx and GEO, drugs proposed should be considered as potential treatments for aging.","summary":"auto","gpt_summary":""}}
  },
  {
    id: '9cb775b8-bbee-3062-e2db-8b04125ef3d3',
    label: 'Use Case 4: Identify the Tissue Activity for a TF based on its Targets',
    description: `Given a Transcription Factor, we collect its targets by various resources and then Enrich the set of consensus targets against GTEx Tissue expression. The result is ranked tissues potentially regulated by the transcription factor.`,
    gpt_summary: 'To start the workflow, KLF4 is selected as the search term. Matching terms from three libraries are identified: ENCODE TF ChIP-seq 2015, ChEA 2022, and ARCHS4 TF Co-Expression. The matching terms are assembled into gene sets, which are then combined into one library. A consensus gene set is created by retaining genes that appear in at least two sets. The gene set is enriched against the GTEx Tissues V8 2023 library to identify statistically significant GTEx Tissue Signatures. The libraries used in this workflow are: \n1. An integrated encyclopedia of DNA elements in the human genome. Nature vol. 489 57–74 (2012). doi:10.1038/nature11247 \n2. Keenan, A. B. et al. ChEA3: transcription factor enrichment analysis by orthogonal omics integration. Nucleic Acids Research vol. 47 W212–W224 (2019). doi:10.1093/nar/gkz446 \n3. Lachmann, A. et al. Massive mining of publicly available RNA-seq data from human and mouse. Nature Communications vol. 9 (2018). doi:10.1038/s41467-018-03751-6 \n4. Lonsdale, J. et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics vol. 45 580–585 (2013). doi:10.1038/ng.2653',
    published: 'Mar 16, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CC-BY-NC-SA-4.0',
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
    workflow: {"data":{"7a70761b-5310-ae85-7095-21adf12e5520":{"type":"Term[Gene]","value":"KLF4"}},"workflow":[{"id":"3a2b64f6-8657-d7ef-2fc6-d7e2f8f74708","type":"Input[Gene]","data":{"id":"7a70761b-5310-ae85-7095-21adf12e5520"}},{"id":"9990ab27-dd9f-4511-47ab-fc4113727949","type":"EnrichrTermSearch[Gene]","inputs":{"term":{"id":"3a2b64f6-8657-d7ef-2fc6-d7e2f8f74708"}}},{"id":"20ec4cb4-2bac-62b8-1528-c037f6047782","type":"ExtractEnrichrTermSearch[ENCODE_TF_ChIP-seq_2015]","inputs":{"searchResults":{"id":"9990ab27-dd9f-4511-47ab-fc4113727949"}}},{"id":"5477675d-88dc-75ee-5b03-5daea612536f","type":"EnrichrSetTToGMT[Gene]","inputs":{"enrichrset":{"id":"20ec4cb4-2bac-62b8-1528-c037f6047782"}}},{"id":"681a371d-57d6-c256-b0aa-b12fa8ade2c2","type":"ExtractEnrichrTermSearch[ChEA_2022]","inputs":{"searchResults":{"id":"9990ab27-dd9f-4511-47ab-fc4113727949"}}},{"id":"e055ec8d-319e-ee8e-c375-f6e05131aaf7","type":"EnrichrSetTToGMT[Gene]","inputs":{"enrichrset":{"id":"681a371d-57d6-c256-b0aa-b12fa8ade2c2"}}},{"id":"ea13ffde-0b13-57ed-b61b-20ca1a05bd60","type":"ExtractEnrichrTermSearch[ARCHS4_TFs_Coexp]","inputs":{"searchResults":{"id":"9990ab27-dd9f-4511-47ab-fc4113727949"}}},{"id":"a7047af0-6aa0-5161-084c-0d35d72d5f87","type":"EnrichrSetTToGMT[Gene]","inputs":{"enrichrset":{"id":"ea13ffde-0b13-57ed-b61b-20ca1a05bd60"}}},{"id":"3d70ee2b-f8ec-c21b-a9b3-155b5421b9ce","type":"GMTConcatenate","inputs":{"gmts:0":{"id":"5477675d-88dc-75ee-5b03-5daea612536f"},"gmts:1":{"id":"e055ec8d-319e-ee8e-c375-f6e05131aaf7"},"gmts:2":{"id":"a7047af0-6aa0-5161-084c-0d35d72d5f87"}}},{"id":"37a07fba-bf8a-eca8-35cf-84362ec4b04e","type":"GMTConsensus","inputs":{"gmt":{"id":"3d70ee2b-f8ec-c21b-a9b3-155b5421b9ce"}}},{"id":"fd13871c-66c6-58be-7796-fb4b7212528c","type":"EnrichrGenesetSearch","inputs":{"geneset":{"id":"37a07fba-bf8a-eca8-35cf-84362ec4b04e"}}},{"id":"44e3f8d4-003f-2937-ab32-c6a393626ea5","type":"ExtractEnrichrGenesetSearch[GTEx_Tissues_V8_2023]","inputs":{"searchResults":{"id":"fd13871c-66c6-58be-7796-fb4b7212528c"}}}],"metadata":{"id":"f30d470a-597c-9f84-f002-1cd06f7d2bf3","title":"Use Case 4: Identify the Tissue Activity for a TF based on its Targets","description":"Given a Transcription Factor, we collect its targets by various resources and then Enrich the set of consensus targets against GTEx Tissue expression. The result is ranked tissues potentially regulated by the transcription factor.","summary":"auto","gpt_summary":""}},
  },
  {
    id: '77805f8b-978f-2ab0-3268-d40c8a06e692',
    label: 'Use Case 5: Small Molecules to Induce a Biological Process',
    description: `We identify genes associated with a biological process from human, mouse phenotypes, KEGG pathways and GO gene set libraries. We then find Consensus LINCS compounds which upregulate these genes resulting in a ranked listing of drug candidates for inducing the biological process.`,
    gpt_summary: `To start the workflow, the search term "Autophagy" was selected. Gene sets with labels containing Autophagy were obtained from Enrichr [1]. These matching terms from the MGI Mammalian Phenotype Level 4 2019 library [2] were then collected to create a collection of gene sets. A GMT file was extracted from the Enrichr results for MGI_Mammalian_Phenotype_Level_4_2019. All the identified gene sets were combined using the union set operation.

Next, reversers and mimickers from over 1 million signatures were identified using SigCom LINCS [3]. The resolved drugs from the LINCS L1000 Chemical Perturbagens library were also included. Matching terms from the KEGG 2021 Human library [4] were gathered to create another collection of gene sets. A GMT file was extracted from the Enrichr results for KEGG_2021_Human. Once again, the identified gene sets were combined using the union set operation.

Reversers and mimickers from over 1 million signatures were identified using SigCom LINCS [3]. Additionally, matching terms from the GO Biological Process 2021 library [5] were collected to create a collection of gene sets. A GMT file was extracted from the Enrichr results for GO_Biological_Process_2021. All the identified gene sets were combined using the union set operation.

Once again, reversers and mimickers from over 1 million signatures were identified using SigCom LINCS [3]. The resolved drugs from the LINCS L1000 Chemical Perturbagens library were included. Mean scores were computed. The drugs were filtered by FDA Approved Drugs with the help of PubChem APIs [6].

References:
1. Xie, Z. et al. Gene Set Knowledge Discovery with Enrichr. Current Protocols vol. 1 (2021). doi:10.1002/cpz1.90
2. Blake, J. A. et al. Mouse Genome Database (MGD): Knowledgebase for mouse–human comparative biology. Nucleic Acids Research vol. 49 D981–D987 (2020). doi:10.1093/nar/gkaa1083
3. Evangelista, J. E. et al. SigCom LINCS: data and metadata search engine for a million gene expression signatures. Nucleic Acids Research vol. 50 W697–W709 (2022). doi:10.1093/nar/gkac328
4. Kanehisa, M., Furumichi, M., Sato, Y., Kawashima, M. & Ishiguro-Watanabe, M. KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Research vol. 51 D587–D592 (2022). doi:10.1093/nar/gkac963
5. Ashburner, M. et al. Gene Ontology: tool for the unification of biology. Nature Genetics vol. 25 25–29 (2000). doi:10.1038/75556
6. Kim, S. et al. PubChem 2023 update. Nucleic Acids Research vol. 51 D1373–D1380 (2022). doi:10.1093/nar/gkac956`,
    published: 'July 26, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CC-BY-NC-SA-4.0',
    dataSources: [
      'HPO',
      'KOMP',
      'KEGG',
      'WikiPathways',
      'GO',
      'LINCS',
    ],
    inputs: [
      PhenotypeTerm,
    ],
    outputs: [
      ScoredDrugs,
    ],
    workflow: {"data":{"ac647b33-1312-0865-c3b1-0ec30840c28f":{"type":"Term[Pathway]","value":"Autophagy"}},"workflow":[{"id":"7fd69f96-85dc-68cb-61c9-077255221f05","type":"Input[Pathway]","data":{"id":"ac647b33-1312-0865-c3b1-0ec30840c28f"}},{"id":"56f988e3-d6be-b2b6-daed-10c8f791a616","type":"EnrichrTermSearch[Pathway]","inputs":{"term":{"id":"7fd69f96-85dc-68cb-61c9-077255221f05"}}},{"id":"72148de6-50ab-1ce6-5c8a-5f40b1febfd0","type":"ExtractEnrichrTermSearch[MGI_Mammalian_Phenotype_Level_4_2019]","inputs":{"searchResults":{"id":"56f988e3-d6be-b2b6-daed-10c8f791a616"}}},{"id":"0fb63e8e-70bb-c841-b01b-68ce535375c1","type":"EnrichrSetTToGMT[Phenotype]","inputs":{"enrichrset":{"id":"72148de6-50ab-1ce6-5c8a-5f40b1febfd0"}}},{"id":"6e708b71-ba56-1382-cef3-e15865e9eb53","type":"GMTUnion","inputs":{"gmt":{"id":"0fb63e8e-70bb-c841-b01b-68ce535375c1"}}},{"id":"4915a324-2704-4e78-ee50-b41824943f46","type":"SigComLINCSGeneSetSearch","inputs":{"genes":{"id":"6e708b71-ba56-1382-cef3-e15865e9eb53"}}},{"id":"85df101b-fe48-5145-0b60-76c0e4afe7cd","type":"ExtractSigComLINCSGeneSetSearch[LINCS L1000 Chemical Perturbagens]","inputs":{"searchResults":{"id":"4915a324-2704-4e78-ee50-b41824943f46"}}},{"id":"17e0f3c8-778b-e4b6-08ac-3f8db370b87e","type":"ExtractEnrichrTermSearch[KEGG_2021_Human]","inputs":{"searchResults":{"id":"56f988e3-d6be-b2b6-daed-10c8f791a616"}}},{"id":"86789631-1c23-0fd5-b991-abde941d5fcc","type":"EnrichrSetTToGMT[Pathway]","inputs":{"enrichrset":{"id":"17e0f3c8-778b-e4b6-08ac-3f8db370b87e"}}},{"id":"62ccf69a-3194-da10-7f2e-76ab3848f731","type":"GMTUnion","inputs":{"gmt":{"id":"86789631-1c23-0fd5-b991-abde941d5fcc"}}},{"id":"112be381-243e-bdee-28d1-5e69fd71e553","type":"SigComLINCSGeneSetSearch","inputs":{"genes":{"id":"62ccf69a-3194-da10-7f2e-76ab3848f731"}}},{"id":"bae2fe50-250e-f712-c33e-a331faaf7a8d","type":"ExtractEnrichrTermSearch[GO_Biological_Process_2021]","inputs":{"searchResults":{"id":"56f988e3-d6be-b2b6-daed-10c8f791a616"}}},{"id":"7d70271f-ee6d-4d90-5003-f680455308c2","type":"EnrichrSetTToGMT[Pathway]","inputs":{"enrichrset":{"id":"bae2fe50-250e-f712-c33e-a331faaf7a8d"}}},{"id":"7cb48105-948a-af24-452f-76abd3f60205","type":"GMTUnion","inputs":{"gmt":{"id":"7d70271f-ee6d-4d90-5003-f680455308c2"}}},{"id":"c3229ad6-880f-c66a-b493-2ff41622c7da","type":"SigComLINCSGeneSetSearch","inputs":{"genes":{"id":"7cb48105-948a-af24-452f-76abd3f60205"}}},{"id":"4e5f6cf3-8175-4800-dc9d-fd021a7a9cf4","type":"ExtractSigComLINCSGeneSetSearch[LINCS L1000 Chemical Perturbagens]","inputs":{"searchResults":{"id":"112be381-243e-bdee-28d1-5e69fd71e553"}}},{"id":"d7bc7f03-c05d-6d7d-56b5-37bce6a18d6f","type":"ExtractSigComLINCSGeneSetSearch[LINCS L1000 Chemical Perturbagens]","inputs":{"searchResults":{"id":"c3229ad6-880f-c66a-b493-2ff41622c7da"}}},{"id":"023f83a8-a124-e755-f7c4-18273fea5271","type":"MeanScoredT[Drug]","inputs":{"scored:0":{"id":"85df101b-fe48-5145-0b60-76c0e4afe7cd"},"scored:1":{"id":"4e5f6cf3-8175-4800-dc9d-fd021a7a9cf4"},"scored:2":{"id":"d7bc7f03-c05d-6d7d-56b5-37bce6a18d6f"}}},{"id":"f0bd5771-7533-5e3a-61ce-28f259b09186","type":"FilterFDAApprovedDrugs[Scored[Drug]]","inputs":{"drugs":{"id":"023f83a8-a124-e755-f7c4-18273fea5271"}}}],"metadata":{"id":"492efde0-c90e-14a8-7859-e3d618f6f261","title":"Use Case 5: Small Molecules to Induce a Biological Process","description":"We identify genes associated with a biological process from human, mouse phenotypes, KEGG pathways and GO gene set libraries. We then find Consensus LINCS compounds which upregulate these genes resulting in a ranked listing of drug candidates for inducing the biological process.","summary":"auto","gpt_summary":""}},
  },
  {
    id: '78db3b40-bf62-4415-ac37-032e8c089859',
    label: 'Use Case 6: CFDE Knowledge about a Variant',
    description: `We query several CFDE data sources for information about the variant provided and about the closest gene to that variant.`,
    gpt_summary: `To begin the workflow, the search term chr10:g.3823823G>A is selected. The closest gene to the variant is determined using MyVariant.info [1]. Next, RNA-seq-like LINCS L1000 Signatures [2] are used to visualize the expression of KLF6, either mimicking or reversing it. The median expression of KLF6 is obtained from the GTEx Portal [3] using the portal's API. Finally, a bar plot is created to visualize the level of expression across scored tissues.

References:
1. Lelong, S. et al. BioThings SDK: a toolkit for building high-performance data APIs in biomedical research. Bioinformatics vol. 38 2077–2079 (2022). doi:10.1093/bioinformatics/btac017
2. Evangelista, J. E. et al. SigCom LINCS: data and metadata search engine for a million gene expression signatures. Nucleic Acids Research vol. 50 W697–W709 (2022). doi:10.1093/nar/gkac328
3. Lonsdale, J. et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics vol. 45 580–585 (2013). doi:10.1038/ng.2653`,
    published: 'July 26, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CC-BY-NC-SA-4.0',
    dataSources: [
      'BioThings',
      'LINCS',
      'GTEx',
    ],
    inputs: [
      VariantTerm,
    ],
    outputs: [
      LINCSL1000ReverseSearchDashboard,
      PlotlyPlot,
    ],
    workflow: {"data":{"a8aa5827-d86b-c54e-31aa-eeed8c4387a7":{"type":"Term[Variant]","value":"chr10:g.3823823G>A"}},"workflow":[{"id":"10532149-325b-7997-53ff-2a423345c73a","type":"Input[Variant]","data":{"id":"a8aa5827-d86b-c54e-31aa-eeed8c4387a7"}},{"id":"f959b7e7-b3c1-f966-e125-0478bf0855b4","type":"GeneTermFromVariantTerm","inputs":{"variant":{"id":"10532149-325b-7997-53ff-2a423345c73a"}}},{"id":"eeab802d-9546-8fe7-a9e7-09cfd7773d2b","type":"LINCSL1000ReverseSearch","inputs":{"gene":{"id":"f959b7e7-b3c1-f966-e125-0478bf0855b4"}}},{"id":"2c43be05-c759-872d-30c1-fb42341afdfd","type":"GTExTissueExpressionFromGene","inputs":{"gene":{"id":"f959b7e7-b3c1-f966-e125-0478bf0855b4"}}},{"id":"7993426b-970b-f949-4e0c-a9d4aa4883f4","type":"BarplotFrom[Scored[Tissue]]","inputs":{"terms":{"id":"2c43be05-c759-872d-30c1-fb42341afdfd"}}}],"metadata":{"id":"9040d19b-3e4a-fd68-67c9-46aa1b391dda","title":"Use Case 6: CFDE Knowledge about a Variant","description":"We query several CFDE data sources for information about the variant provided and about the closest gene to that variant.","summary":"auto","gpt_summary":""}},
  },
  {
    id: '60086ff7-2a49-195a-9682-146e763e32cc',
    label: 'Use Case 6: CFDE Knowledge about a Gene',
    description: `We query several CFDE data sources for information about the gene provided.`,
    gpt_summary: `To start the workflow, the search term selected is KLF6. Next, we visualize RNA-seq-like LINCS L1000 Signatures [1] that either mimic or reverse the expression of KLF6. To obtain the median expression of KLF6, we utilize the GTEx Portal [2] and its API. To better understand the expression levels across different tissues, a bar plot is created.

1. Evangelista, J. E. et al. SigCom LINCS: data and metadata search engine for a million gene expression signatures. Nucleic Acids Research vol. 50 W697–W709 (2022). doi:10.1093/nar/gkac328
2. Lonsdale, J. et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics vol. 45 580–585 (2013). doi:10.1038/ng.2653`,
    published: 'July 26, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CC-BY-NC-SA-4.0',
    dataSources: [
      'BioThings',
      'LINCS',
      'GTEx',
    ],
    inputs: [
      GeneTerm,
    ],
    outputs: [
      LINCSL1000ReverseSearchDashboard,
      PlotlyPlot,
    ],
    workflow: {"data":{"caf3189c-df8d-bd9d-38f1-cdbe6091fb78":{"type":"Term[Gene]","value":"KLF6"}},"workflow":[{"id":"dbc3db82-3dbe-facf-f88b-bff608fd99b6","type":"Input[Gene]","data":{"id":"caf3189c-df8d-bd9d-38f1-cdbe6091fb78"}},{"id":"9ae1b029-8092-11e4-b8e2-3987a6ba4047","type":"LINCSL1000ReverseSearch","inputs":{"gene":{"id":"dbc3db82-3dbe-facf-f88b-bff608fd99b6"}}},{"id":"cb94ed51-6c8f-cdb5-2953-8cda38984696","type":"GTExTissueExpressionFromGene","inputs":{"gene":{"id":"dbc3db82-3dbe-facf-f88b-bff608fd99b6"}}},{"id":"cea6af28-19ed-0a60-05a7-8d9830d01559","type":"BarplotFrom[Scored[Tissue]]","inputs":{"terms":{"id":"cb94ed51-6c8f-cdb5-2953-8cda38984696"}}}],"metadata":{"id":"87d55f15-3b5c-ac33-779f-2348ecbfd24f","title":"Use Case 6: CFDE Knowledge about a Gene","description":"We query several CFDE data sources for information about the gene provided.","summary":"auto","gpt_summary":""}},
  },
  {
    id: '7103143b-f178-7aa7-7753-9d4652ddc931',
    label: 'Use Case 7: Variant Expression in Tumor/Healthy',
    description: `We construct a joint plot showing how the variant's closest gene is expressed in tumors from KidsFirst and healthy human tissue from GTEx.`,
    gpt_summary: `The workflow begins by selecting the search term chr2:g.39417578C>G. The closest gene to the variant was found using MyVariant.info [1]. CDKL4 gene expression in tumors was queried from the Open Pediatric Cancer Atlas API [2]. The median expression of CDKL4 was obtained from the GTEx Portal [3] using the portal's API. To visualize the level of expression across tumor gene expression, a bar plot was created.

References:
1. Lelong, S. et al. BioThings SDK: A toolkit for building high-performance data APIs in biomedical research. Bioinformatics, vol. 38, 2077-2079 (2022). doi:10.1093/bioinformatics/btac017
2. Shapiro, J. A. et al. OpenPBTA: The Open Pediatric Brain Tumor Atlas. Cell Genomics 100340 (2023). doi:10.1016/j.xgen.2023.100340
3. Lonsdale, J. et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics, vol. 45, 580-585 (2013). doi:10.1038/ng.2653`,
    published: 'July 26, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CC-BY-NC-SA-4.0',
    dataSources: [
      'BioThings',
      'KidsFirst',
      'GTEx',
    ],
    inputs: [
      VariantTerm,
    ],
    outputs: [
      PlotlyPlot,
    ],
    workflow: {"data":{"8c6cd90e-b7c4-e411-8d87-96f8aed9c7a9":{"type":"Term[Variant]","value":"chr2:g.39417578C>G"}},"workflow":[{"id":"fb0d7725-ce39-b91c-a346-3a2f780e3c6b","type":"Input[Variant]","data":{"id":"8c6cd90e-b7c4-e411-8d87-96f8aed9c7a9"}},{"id":"6a2e136d-80c0-8601-a934-4069e7732e8f","type":"GeneTermFromVariantTerm","inputs":{"variant":{"id":"fb0d7725-ce39-b91c-a346-3a2f780e3c6b"}}},{"id":"246abb4b-3bea-6f40-fdf0-328c8c780745","type":"KFTumorExpressionFromGene","inputs":{"gene":{"id":"6a2e136d-80c0-8601-a934-4069e7732e8f"}}},{"id":"3ffa28f8-9622-d550-b512-2524fd26ea85","type":"GTExTissueExpressionFromGene","inputs":{"gene":{"id":"6a2e136d-80c0-8601-a934-4069e7732e8f"}}},{"id":"761a3dc2-b523-ebff-d109-02cd2c1b3669","type":"BarplotFrom[[TumorGeneExpression]]","inputs":{"other_terms":{"id":"3ffa28f8-9622-d550-b512-2524fd26ea85"},"terms":{"id":"246abb4b-3bea-6f40-fdf0-328c8c780745"}}}],"metadata":{"id":"29816ef1-3494-a461-73a4-f5a0410a8320","title":"Use Case 7: Variant Expression in Tumor/Healthy","description":"We construct a joint plot showing how the variant's closest gene is expressed in tumors from KidsFirst and healthy human tissue from GTEx.","summary":"auto","gpt_summary":""}},
  },
//   {
//     id: '703a150e-ee10-2a26-418a-89a69e1a82e7',
//     label: 'Use Case 8: Associations between 2 Variants',
//     description: `Given two variants, we find their closest genes and present combined knowledge about them including their expression in tumors & healthy tissue, interactions between the two proteins, and gene sets containing the two genes.`,
//     gpt_summary: `The workflow begins by selecting chr10:g.3823823G>A as the search term. Following that, chr2:g.39417578C>G is chosen as the search term. The closest gene to each variant is determined using MyVariant.info [1]. The specified genes are then combined into a gene set. This gene set is submitted to Enrichr [2]. Gene expression data for CDKL4 in tumors is obtained from the Open Pediatric Cancer Atlas API [3]. Median expression of CDKL4 is also obtained from the GTEx Portal [4] using the portal's API. A bar plot is created to visualize the level of expression across tumor gene expression for CDKL4. Similarly, median expression of KLF6 is obtained from the GTEx Portal [4] using the portal's API. Gene expression data for KLF6 in tumors is queried from the Open Pediatric Cancer Atlas API [3]. Another bar plot is created to visualize the level of expression across tumor gene expression for KLF6.

// References:
// 1. Lelong, S. et al. BioThings SDK: a toolkit for building high-performance data APIs in biomedical research. Bioinformatics vol. 38 2077–2079 (2022). doi:10.1093/bioinformatics/btac017
// 2. Xie, Z. et al. Gene Set Knowledge Discovery with Enrichr. Current Protocols vol. 1 (2021). doi:10.1002/cpz1.90
// 3. Shapiro, J. A. et al. OpenPBTA: The Open Pediatric Brain Tumor Atlas. Cell Genomics 100340 (2023) doi:10.1016/j.xgen.2023.100340. doi:10.1016/j.xgen.2023.100340
// 4. Lonsdale, J. et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics vol. 45 580–585 (2013). doi:10.1038/ng.2653`,
//     published: 'July 26, 2023',
//     version: '1.0.0',
//     authors: ['CFDE Playbook Partnership'],
//     licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
//     license: 'CC-BY-NC-SA-4.0',
//     dataSources: [
//       'KidsFirst',
//       'GTEx',
//       'Enrichr',
//       'BioThings',
//     ],
//     inputs: [
//       VariantTerm,
//     ],
//     outputs: [
//       PlotlyPlot,
//     ],
//     workflow: {"data":{"a8aa5827-d86b-c54e-31aa-eeed8c4387a7":{"type":"Term[Variant]","value":"chr10:g.3823823G>A"},"8c6cd90e-b7c4-e411-8d87-96f8aed9c7a9":{"type":"Term[Variant]","value":"chr2:g.39417578C>G"}},"workflow":[{"id":"10532149-325b-7997-53ff-2a423345c73a","type":"Input[Variant]","data":{"id":"a8aa5827-d86b-c54e-31aa-eeed8c4387a7"}},{"id":"fb0d7725-ce39-b91c-a346-3a2f780e3c6b","type":"Input[Variant]","data":{"id":"8c6cd90e-b7c4-e411-8d87-96f8aed9c7a9"}},{"id":"f959b7e7-b3c1-f966-e125-0478bf0855b4","type":"GeneTermFromVariantTerm","inputs":{"variant":{"id":"10532149-325b-7997-53ff-2a423345c73a"}}},{"id":"6a2e136d-80c0-8601-a934-4069e7732e8f","type":"GeneTermFromVariantTerm","inputs":{"variant":{"id":"fb0d7725-ce39-b91c-a346-3a2f780e3c6b"}}},{"id":"8b4b2a6b-6ab6-56c2-0f34-cd9e08b66988","type":"TermToSet[Term[Gene]]","inputs":{"terms:0":{"id":"f959b7e7-b3c1-f966-e125-0478bf0855b4"},"terms:1":{"id":"6a2e136d-80c0-8601-a934-4069e7732e8f"}}},{"id":"ff97e368-5587-f981-e92c-585c38d3a324","type":"EnrichrGenesetSearch","inputs":{"geneset":{"id":"8b4b2a6b-6ab6-56c2-0f34-cd9e08b66988"}}},{"id":"246abb4b-3bea-6f40-fdf0-328c8c780745","type":"KFTumorExpressionFromGene","inputs":{"gene":{"id":"6a2e136d-80c0-8601-a934-4069e7732e8f"}}},{"id":"3ffa28f8-9622-d550-b512-2524fd26ea85","type":"GTExTissueExpressionFromGene","inputs":{"gene":{"id":"6a2e136d-80c0-8601-a934-4069e7732e8f"}}},{"id":"761a3dc2-b523-ebff-d109-02cd2c1b3669","type":"BarplotFrom[[TumorGeneExpression]]","inputs":{"other_terms":{"id":"3ffa28f8-9622-d550-b512-2524fd26ea85"},"terms":{"id":"246abb4b-3bea-6f40-fdf0-328c8c780745"}}},{"id":"2c43be05-c759-872d-30c1-fb42341afdfd","type":"GTExTissueExpressionFromGene","inputs":{"gene":{"id":"f959b7e7-b3c1-f966-e125-0478bf0855b4"}}},{"id":"f1cdfc2e-83c9-07f5-09c8-c9ab5690a800","type":"KFTumorExpressionFromGene","inputs":{"gene":{"id":"f959b7e7-b3c1-f966-e125-0478bf0855b4"}}},{"id":"181241ce-d1e9-46bf-5cb3-4b4c6b6d3c51","type":"BarplotFrom[[TumorGeneExpression]]","inputs":{"other_terms":{"id":"2c43be05-c759-872d-30c1-fb42341afdfd"},"terms":{"id":"f1cdfc2e-83c9-07f5-09c8-c9ab5690a800"}}}],"metadata":{"id":"c6eadeb9-7456-6913-ed75-bed0a72c1d65","title":"Use Case 8: Associations between 2 Variants","description":"","summary":"auto","gpt_summary":""}},
//   },
//   {
//     id: 'ec5a2d97-fca2-df40-578a-9c84129963bf',
//     label: 'Use Case 8: Associations between 2 Genes',
//     description: `We present combined knowledge about the two given genes including their expression in tumors & healthy tissue, interactions between the two proteins, and gene sets containing the two genes.`,
//     gpt_summary: `The workflow begins by selecting ACE2 and STAT3 as the search terms. The specified genes are then combined into a gene set. This gene set is submitted to Enrichr [1]. The median expression of ACE2 is obtained from the GTEx Portal [2] using the portal's API. The gene expression of ACE2 in tumors is queried from the Open Pediatric Cancer Atlas API [3]. A bar plot is created to visualize the level of expression across tumor gene expression for ACE2. Similarly, the median expression of STAT3 is obtained from the GTEx Portal [2] using the portal's API. The gene expression of STAT3 in tumors is queried from the Open Pediatric Cancer Atlas API [3]. Again, a bar plot is created to visualize the level of expression across tumor gene expression for STAT3.

//     References:
//     1. Xie, Z. et al. Gene Set Knowledge Discovery with Enrichr. Current Protocols vol. 1 (2021). doi:10.1002/cpz1.90
//     2. Lonsdale, J. et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics vol. 45 580–585 (2013). doi:10.1038/ng.2653
//     3. Shapiro, J. A. et al. OpenPBTA: The Open Pediatric Brain Tumor Atlas. Cell Genomics 100340 (2023) doi:10.1016/j.xgen.2023.100340. doi:10.1016/j.xgen.2023.100340.`,
//     published: 'July 26, 2023',
//     version: '1.0.0',
//     authors: ['CFDE Playbook Partnership'],
//     licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
//     license: 'CC-BY-NC-SA-4.0',
//     dataSources: [
//       'KidsFirst',
//       'GTEx',
//       'Enrichr',
//     ],
//     inputs: [
//       VariantTerm,
//     ],
//     outputs: [
//       PlotlyPlot,
//     ],
//     workflow: {"data":{"17165506-14be-2602-ba46-bc5f96d1ae6d":{"type":"Term[Gene]","value":"ACE2"},"4d8b4e54-b938-fcdd-d6d2-1fae7a0e627c":{"type":"Term[Gene]","value":"STAT3"}},"workflow":[{"id":"85f38624-d0be-9c74-c9d8-6cf3a903f53c","type":"Input[Gene]","data":{"id":"17165506-14be-2602-ba46-bc5f96d1ae6d"}},{"id":"1d6aaca9-1aae-d897-e954-c22ce89671e8","type":"Input[Gene]","data":{"id":"4d8b4e54-b938-fcdd-d6d2-1fae7a0e627c"}},{"id":"7e38312e-6c0e-8cae-f264-8738774c2d42","type":"TermToSet[Term[Gene]]","inputs":{"terms:0":{"id":"85f38624-d0be-9c74-c9d8-6cf3a903f53c"},"terms:1":{"id":"1d6aaca9-1aae-d897-e954-c22ce89671e8"}}},{"id":"94c9f891-c6e7-e95e-b906-ec96f340da66","type":"EnrichrGenesetSearch","inputs":{"geneset":{"id":"7e38312e-6c0e-8cae-f264-8738774c2d42"}}},{"id":"681c1040-cd86-6acc-1c96-ea82c54bde0e","type":"GTExTissueExpressionFromGene","inputs":{"gene":{"id":"85f38624-d0be-9c74-c9d8-6cf3a903f53c"}}},{"id":"eba1261b-b32d-25af-d7de-78c4cf5ec455","type":"KFTumorExpressionFromGene","inputs":{"gene":{"id":"85f38624-d0be-9c74-c9d8-6cf3a903f53c"}}},{"id":"a0d2698d-123e-6e42-306b-51421e1c87c4","type":"BarplotFrom[[TumorGeneExpression]]","inputs":{"other_terms":{"id":"681c1040-cd86-6acc-1c96-ea82c54bde0e"},"terms":{"id":"eba1261b-b32d-25af-d7de-78c4cf5ec455"}}},{"id":"392bef46-2b9c-d4e8-5b4e-ef63099bff88","type":"GTExTissueExpressionFromGene","inputs":{"gene":{"id":"1d6aaca9-1aae-d897-e954-c22ce89671e8"}}},{"id":"4f9e9fe9-93fb-c742-b748-913d86c490a4","type":"KFTumorExpressionFromGene","inputs":{"gene":{"id":"1d6aaca9-1aae-d897-e954-c22ce89671e8"}}},{"id":"e918d85c-f942-0758-00d2-fdd7c65001d0","type":"BarplotFrom[[TumorGeneExpression]]","inputs":{"other_terms":{"id":"392bef46-2b9c-d4e8-5b4e-ef63099bff88"},"terms":{"id":"4f9e9fe9-93fb-c742-b748-913d86c490a4"}}}],"metadata":{"id":"e60c11d9-ff31-af9e-a790-cc5a9df8f064","title":"Use Case 8: Associations between 2 Genes","description":"","summary":"auto","gpt_summary":""}},
//   },
  // {
  //   id: 'cc338d00-d22e-08b1-3b74-21710dec3847',
  //   label: 'Use Case 9: Identifying regulatory relationships between genes, regulatory regions, and variants',
  //   description: `Using exRNA APIs, we identify tissue and allele specific regulatory information to perform regulatory variant burden testing.`,
  //   gpt_summary: `The workflow begins by selecting CA13203640 as the search term. Additional information about the variant CA13203640 is then resolved. Regulatory elements associated with the variant CA13203640 are also resolved. Further information about the regulatory element EH38E2924876 is obtained. Genes that are linked to the regulatory element EH38E2924876 are determined. Variants that are linked to the regulatory element EH38E2924876 are also resolved. Allele specific evidences for the variant CA13203640 are then resolved. Lastly, xQTL evidence data for the variant CA13203640 is obtained.`,
  //   published: 'July 26, 2023',
  //   version: '1.0.0',
  //   authors: ['CFDE Playbook Partnership'],
  //   licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
  //   license: 'CC-BY-NC-SA-4.0',
  //   dataSources: [
  //     'ClinGen',
  //     'STRING',
  //     'KEGG',
  //   ],
  //   inputs: [
  //     VariantTerm,
  //   ],
  //   outputs: [
  //     RegulatoryElementSet,
  //   ],
  //   workflow: {"data":{"2769d61e-6fed-2989-060f-7c3dd935b2fb":{"type":"Term[Variant]","value":"CA13203640"}},"workflow":[{"id":"35d0a1ed-1bb5-fd4d-a351-50095e8c516b","type":"Input[Variant]","data":{"id":"2769d61e-6fed-2989-060f-7c3dd935b2fb"}},{"id":"0297e99b-c758-9137-95c5-23a9c06f1bfb","type":"VariantInfoFromVariantTerm","inputs":{"variant":{"id":"35d0a1ed-1bb5-fd4d-a351-50095e8c516b"}}},{"id":"b75bde78-8b11-ddf0-386a-d0d9d1fc5ee2","type":"GetRegulatoryElementForThisVariant","inputs":{"variantInfo":{"id":"0297e99b-c758-9137-95c5-23a9c06f1bfb"}}},{"id":"48f040ff-60ae-7553-1803-61313c1b3240","type":"RegElementInfoFromRegElementTerm","inputs":{"regulatoryElement":{"id":"b75bde78-8b11-ddf0-386a-d0d9d1fc5ee2"}}},{"id":"28f2ffac-71ab-f961-9ae9-02ce2a8d09b9","type":"GetGenesForRegulatoryElementInfo","inputs":{"regElemInfo":{"id":"48f040ff-60ae-7553-1803-61313c1b3240"}}},{"id":"2c815b01-388d-e5ba-b4e9-dc43902e7c49","type":"GetVariantListForRegulatoryElementInfo","inputs":{"regElemInfo":{"id":"48f040ff-60ae-7553-1803-61313c1b3240"}}},{"id":"0d40f27a-e532-d46d-5581-6d64b4f31d24","type":"GetAlleleSpecificEvidencesForThisVariant","inputs":{"variantInfo":{"id":"0297e99b-c758-9137-95c5-23a9c06f1bfb"}}},{"id":"3420df8c-8b71-8d07-f3c7-9b61be9b4061","type":"GetxQTL_EvidencesDataForVariantInfo","inputs":{"variantInfo":{"id":"0297e99b-c758-9137-95c5-23a9c06f1bfb"}}}],"metadata":{"id":"4a760fc1-bca1-b5d1-0aee-11d4e762a66a","title":"Use Case 9: Identifying regulatory relationships between genes, regulatory regions, and variants","description":"Using exRNA APIs, we identify tissue and allele specific regulatory information to perform regulatory variant burden testing.","summary":"auto","gpt_summary":""}},
  // },
  {
    id: '9bad4205-70e2-7f8d-b151-8fbb11e58ea9',
    label: 'Use Case 10: Guilt by Association',
    description: `Given a set of genes, connect the dots (CTD) is performed against protein & pathway graphs to obtain a small subset of highly connected genes and those that are guilty by association.`,
    gpt_summary: `The workflow begins by creating a gene set from Example geneset. Then, CTD is performed on the gene set using STRING and KEGG. From the CTD output, a list of Highly Connected Genes is obtained, as well as a list of Guilty By Association Genes. This process is repeated, resulting in another list of Highly Connected Genes and another list of Guilty By Association Genes.`,
    published: 'July 26, 2023',
    version: '1.0.1',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CC-BY-NC-SA-4.0',
    dataSources: [
      'ClinGen',
      'STRING',
    ],
    inputs: [
      GeneSet,
    ],
    outputs: [
      GeneSet,
    ],
    workflow: {"data":{"60e507db-c00d-d5eb-d014-6498e81aed4f":{"type":"Input[Set[Gene]]","value":{"description":"Example gene set","set":["UTP14A","S100A6","SCAND1","RRP12","CIAPIN1","ADH5","MTERF3","SPR","CHMP4A","UFM1","VAT1","HACD3","RFC5","COTL1","NPRL2","TRIB3","PCCB","TLE1","CD58","BACE2","KDM3A","TARBP1","RNH1","CHAC1","MBNL2","VDAC1","TES","OXA1L","NOP56","HAT1","CPNE3","DNMT1","ARHGAP1","VPS28","EIF2S2","BAG3","CDCA4","NPDC1","RPS6KA1","FIS1","SYPL1","SARS","CDC45","CANT1","HERPUD1","SORBS3","MRPS2","TOR1A","TNIP1","SLC25A46","MAL","EPCAM","HDAC6","CAPN1","TNRC6B","PKD1","RRS1","HP","ANO10","CEP170B","IDE","DENND2D","CAMK2B","ZNF358","RPP38","MRPL19","NUCB2","GNAI1","LSR","ADGRE2","PKMYT1","CDK5R1","ABL1","PILRB","AXIN1","FBXL8","MCF2L","DBNDD1","IGHMBP2","WIPF2","WFS1","OGFOD2","MAPK1IP1L","COL11A1","REG3A","SERPINA1","MYCBP2","PIGK","TCAP","CRADD","ELK1","DNAJB2","ZBTB16","DAZAP1","MAPKAPK2","EDRF1","CRIP1","UCP3","AGR2","P4HA2"]}}},"workflow":[{"id":"df5b091e-7c80-f268-9967-dbaffc68d89c","type":"Input[Set[Gene]]","data":{"id":"60e507db-c00d-d5eb-d014-6498e81aed4f"}},{"id":"0b3cae1b-2e8d-4fab-ad1c-42fb483a5e62","type":"GeneSet_CTD_String","inputs":{"geneset":{"id":"df5b091e-7c80-f268-9967-dbaffc68d89c"}}},{"id":"b838fc0a-ab33-574c-6967-f9d8699925bc","type":"Highly_Connected_Genes","inputs":{"ctdResponseInfo":{"id":"0b3cae1b-2e8d-4fab-ad1c-42fb483a5e62"}}},{"id":"ec5da81b-6bc2-66b9-713c-c1d7de9ec3c6","type":"Guilty_By_Association_Genes","inputs":{"ctdResponseInfo":{"id":"0b3cae1b-2e8d-4fab-ad1c-42fb483a5e62"}}}],"metadata":{"id":"d02a3634-df43-ee76-a23f-7f23b0a16138","title":"Use Case 10: Guilt by Association","description":"Given a set of genes, connect the dots (CTD) is performed against protein & pathway graphs to obtain a small subset of highly connected genes and those that are guilty by association.","summary":"auto","gpt_summary":""}},
  },
  {
    id: '24f9b280-82fa-4fb3-94e9-1ad80c7b9f3e',
    label: 'Use Case 11: Related Proteins/Metabolites across DCCs',
    description: `Using directly related proteins to a gene of interest, we return a slew of related information to the gene protein network from different DCCs.`,
    gpt_summary: `To begin the workflow, the search term RPE is selected. Using the API from StringDB PPI, the gene ID (SYMBOL) is extracted [1]. This extraction generates a list of nodes (Gene Set) from the given StringDB PPI. From over 1 million signatures, reversers and mimickers are identified using SigCom LINCS [2]. The gene set is then submitted to Enrichr [3]. Next, the gene set is searched in the Metabolomics Workbench to identify relevant reactions [4]. Additionally, the gene set is searched in the Metabolomics Workbench to identify associated metabolites [REF]. Finally, the gene set is searched in the Metabolomics Workbench to identify relevant studies related to the genes [5].

Here are the references for the sources mentioned:
1. STRING api, https://string-db.org/cgi/help.pl?subpage=api%23getting-all-the-string-interaction-partners-of-the-protein-set
2. Evangelista, J. E. et al. SigCom LINCS: data and metadata search engine for a million gene expression signatures. Nucleic Acids Research vol. 50 W697–W709 (2022). doi:10.1093/nar/gkac328
3. Xie, Z. et al. Gene Set Knowledge Discovery with Enrichr. Current Protocols vol. 1 (2021). doi:10.1002/cpz1.90
4. The Metabolomics Workbench, https://www.metabolomicsworkbench.org/`,
    published: 'July 27, 2023',
    version: '1.0.0',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CC-BY-NC-SA-4.0',
    dataSources: [
      'LINCS',
      'STRING',
      'ChEA',
      'GTEx',
      'Enrichr',
      'KEGG',
      'GO',
      'MSigDB',
      'Metabolomics',
    ],
    inputs: [
      GeneTerm,
    ],
    outputs: [
      MetGeneRxnTable,
      MetgeneMetaboliteTable,
      MetGeneStudyTable,
    ],
    workflow: {"data":{},"workflow":[{"id":"c9c93e88-c640-63e8-4d5a-487249b0622c","type":"Input[Gene]"},{"id":"f51f891f-01c4-c430-9f7d-ab685df75202","type":"FetchStringDBPPI","inputs":{"gene":{"id":"c9c93e88-c640-63e8-4d5a-487249b0622c"}}},{"id":"b1cd5c16-7c65-d349-3984-d2de48ced74a","type":"StringDBPPI_to_GraphPlot","inputs":{"data":{"id":"f51f891f-01c4-c430-9f7d-ab685df75202"}}},{"id":"aa472914-7301-c427-9b69-3cd5df954cbd","type":"StringDBPPI_to_GeneSet","inputs":{"data":{"id":"f51f891f-01c4-c430-9f7d-ab685df75202"}}},{"id":"ca2e58cb-3999-4f16-0aaf-d8280072d94f","type":"SigComLINCSGeneSetSearch","inputs":{"genes":{"id":"aa472914-7301-c427-9b69-3cd5df954cbd"}}},{"id":"05f8f08d-a91a-c384-23cc-68a55d5f02fc","type":"EnrichrGenesetSearch","inputs":{"geneset":{"id":"aa472914-7301-c427-9b69-3cd5df954cbd"}}},{"id":"bdf62ba1-4685-fb5a-8a0c-8a496d4af069","type":"MetGeneRxnsSet","inputs":{"geneset":{"id":"aa472914-7301-c427-9b69-3cd5df954cbd"}}},{"id":"799e3f4d-7536-44cc-bdbd-6da1f6305963","type":"MetgeneMetabolitesGeneSet","inputs":{"geneset":{"id":"aa472914-7301-c427-9b69-3cd5df954cbd"}}},{"id":"e92bc9f4-7416-8ed1-1c69-55b1caa1460d","type":"MetGeneStudiesGeneSet","inputs":{"geneset":{"id":"aa472914-7301-c427-9b69-3cd5df954cbd"}}}],"metadata":{"id":"ad97b2ae-0b3e-547a-8fdb-fe563a247090","title":"Use Case 11: Related Proteins/Metabolites across DCCs","description":"Using directly related proteins to a gene of interest, we return a slew of related information to the gene protein network from different DCCs.","summary":"auto","gpt_summary":""}},
  },
  // TODO: 12
  {
    id: '63d937db-8082-78f3-4961-57f814a0d0ef',
    label: 'Use Case 13: Novel Cell Surface Targets for Individual Cancer Patients Analyzed with Common Fund Datasets',
    description: `The input to this workflow is a data matrix of gene expression that was collected from a pediatric patient tumor patient from the KidsFirst Common Fund program [1]. The RNA-seq samples are the columns of the matrix, and the rows are the raw expression gene count for all human coding genes (Table 1). This data matrix is fed into TargetRanger [2] to screen for targets which are highly expressed in the tumor but lowly expressed across most healthy human tissues based on gene expression data collected from postmortem patients with RNA-seq by the GTEx Common Fund program [3]. Based on this analysis the gene IMP U3 small nucleolar ribonucleoprotein 3 (IMP3) was selected because it was the top candidate returned from the TargetRanger analysis (Tables 2-3). IMP3 is also commonly called insulin-like growth factor 2 mRNA-binding protein 3 (IGF2BP3). Next, we leverage unique knowledge from various other Common Fund programs to examine various functions and knowledge related to IMP3. First, we queried the LINCS L1000 data [4] from the LINCS program [5] converted into RNA-seq-like LINCS L1000 Signatures [6] using the SigCom LINCS API [7] to identify mimicker or reverser small molecules that maximally impact the expression of IMP3 in human cell lines (Fig. 1, Table 4). In addition, we also queried the LINCS L1000 data to identify single gene CRISPR knockouts that down-regulate the expression of IMP3 (Fig. 1, Table 5). These potential drug targets were filtered using the Common Fund IDG program's list of understudied proteins [8] to produce a set of additional targets (Table 6). Next, IMP3 was searched for knowledge provided by the with the Metabolomics Workbench MetGENE tool [9]. MetGENE aggregates knowledge about pathways, reactions, metabolites, and studies from the Metabolomics Workbench Common Fund supported resource [10]. The Metabolomics Workbench was searched to find associated metabolites linked to IMP3 [10]. Furthermore, we leveraged  the CFDE Linked Data Hub API [11] to list knowledge about regulatory elements associated with IMP3 (Table 6). Finally, the GlyGen database [12] was queried to identify relevant sets of proteins that are the product of the IMP3 genes, as well as known post-translational modifications discovered on IMP3.

References
1. Heath, Allison P., et al. "Gabriella Miller Kids First Data Resource Center: Harmonizing clinical and genomic data to support childhood cancer and structural birth defect research." Cancer Research 79.13_Supplement (2019): 2464-2464.
2. G. B. Marino, et al. GeneRanger and TargetRanger: processed gene and protein expression levels across cells and tissues for target discovery, Nucleic Acids Research. 51, W1, 5 July 2023, Pages W213–W224, https://doi.org/10.1093/nar/gkad399
3. Lonsdale, J. et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics vol. 45 580–585 (2013). doi:10.1038/ng.2653
4. Subramanian A, et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell. 2017 Nov 30;171(6):1437-1452.e17. doi: 10.1016/j.cell.2017.10.049. 
5. Keenan AB, et al. The Library of Integrated Network-Based Cellular Signatures NIH Program: System-Level Cataloging of Human Cells Response to Perturbations. Cell Syst. 2018 Jan 24;6(1):13-24. doi: 10.1016/j.cels.2017.11.001
6. Jeon M, et al. Transforming L1000 profiles to RNA-seq-like profiles with deep learning. BMC Bioinformatics. 2022 Sep 13;23(1):374. doi: 10.1186/s12859-022-04895-5.
7. Evangelista, J. E. et al. SigCom LINCS: data and metadata search engine for a million gene expression signatures. Nucleic Acids Research vol. 50 W697–W709 (2022). doi:10.1093/nar/gkac328
8. Oprea TI, et al. Unexplored therapeutic opportunities in the human genome. Nat Rev Drug Discov. 2018 May;17(5):317-332. doi: 10.1038/nrd.2018.14
9. Srinivasan S, et al. MetGENE: gene-centric metabolomics information retrieval tool. Gigascience. 2022 Dec 28;12:giad089. doi: 10.1093/gigascience/giad089
10. Sud M, et al. Metabolomics Workbench: An international repository for metabolomics data and metadata, metabolite standards, protocols, tutorials and training, and analysis tools. Nucleic Acids Res. 2016 Jan 4;44(D1):D463-70. doi: 10.1093/nar/gkv1042.
11. CFDE Linked Data Hub, https://ldh.genome.network/cfde/ldh/
12. York, W. S. et al. GlyGen: Computational and Informatics Resources for Glycoscience. Glycobiology vol. 30 72–73 (2019). doi:10.1093/glycob/cwz080`,
    gpt_summary: `The process began with the upload of a file, which was then parsed as an input gene count matrix. The next step involved identifying significantly over-expressed genes in comparison to normal tissue in GTEx, resulting in the selection of IMP3 for further investigation.

To identify drugs that down-regulate IMP3 expression, RNA-seq-like LINCS L1000 Chemical Perturbagens were utilized. Additionally, genes that down-regulate IMP3 were identified from RNA-seq-like LINCS L1000 CRISPR Knockouts. The list of genes was then filtered by IDG Understudied Proteins.

The Metabolomics Workbench was used to identify associated metabolites and relevant reactions for IMP3. Regulatory elements were obtained from the CFDE Linked Data Hub. Lastly, the GlyGen database was searched to identify a relevant set of proteins originating from IMP3.`,
    published: 'Mar 30, 2023',
    version: '1.0.2',
    authors: ['CFDE Playbook Partnership'],
    licenseUrl: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
    license: 'CC-BY-NC-SA-4.0',
    dataSources: [
      'KidsFirst',
      'LINCS',
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
      MetGeneRxnTable,
      MetgeneMetaboliteTable,
      VariantSet,
      GlyGenProteinResponseNode,
    ],
    workflow: {"data":{"0a44b5ae-e9ef-a30c-75ed-62fe096374e4":{"type":"GeneCountMatrixFileUpload","value":{"url":"drs://playbook-workflow-builder.cloud/0237e4c2-56c1-49de-a82d-73293c60cc8f","filename":"example.tsv","sha256":"c48cc727487751ecd8b17418c75e9a88aa1ffd4c586d82379737432f00fa5861","size":916155}},"4e00b247-4cf1-02e0-8b77-9acfc5f307ee":{"type":"OneScoredT[Scored[Gene]]","value":"IMP3"}},"workflow":[{"id":"048f537b-ada0-51f3-2ee2-0292cfa4e701","type":"GeneCountMatrixFileUpload","data":{"id":"0a44b5ae-e9ef-a30c-75ed-62fe096374e4"}},{"id":"373e902a-eba7-abaf-c0b0-cf40c9352016","type":"TargetRangerScreenTargets[GTEx_transcriptomics]","inputs":{"input":{"id":"048f537b-ada0-51f3-2ee2-0292cfa4e701"}}},{"id":"8d325af7-e4a3-00fb-16f8-b26a1c97546d","type":"OneScoredT[Scored[Gene]]","inputs":{"scored":{"id":"373e902a-eba7-abaf-c0b0-cf40c9352016"}},"data":{"id":"4e00b247-4cf1-02e0-8b77-9acfc5f307ee"}},{"id":"a74c6ea6-40f2-fa70-70f8-7802c08dba02","type":"LINCSL1000ReverseSearch","inputs":{"gene":{"id":"8d325af7-e4a3-00fb-16f8-b26a1c97546d"}}},{"id":"b91a3ef7-91d2-4e27-e8d0-10c38ecc1543","type":"LINCSL1000ReverseSearchExtract[Drug, Down]","inputs":{"search":{"id":"a74c6ea6-40f2-fa70-70f8-7802c08dba02"}}},{"id":"a87deee7-85ed-2905-dbe0-3170c23c2ca7","type":"LINCSL1000ReverseSearchExtract[Gene, Down]","inputs":{"search":{"id":"a74c6ea6-40f2-fa70-70f8-7802c08dba02"}}},{"id":"f6666b32-d58a-6c22-fee7-a99309ebe58a","type":"IDGFilter[Scored[Gene], all]","inputs":{"input":{"id":"a87deee7-85ed-2905-dbe0-3170c23c2ca7"}}},{"id":"b903ac9e-b77d-d527-5b9c-26e2f58198da","type":"MetGeneSearch","inputs":{"gene":{"id":"8d325af7-e4a3-00fb-16f8-b26a1c97546d"}}},{"id":"2e549c22-53f1-139e-2c1a-8beb62b08607","type":"MetgeneMetabolites","inputs":{"summary":{"id":"b903ac9e-b77d-d527-5b9c-26e2f58198da"}}},{"id":"29b8026f-59ce-7fae-229f-80a0043347a0","type":"MetGeneRxns","inputs":{"summary":{"id":"b903ac9e-b77d-d527-5b9c-26e2f58198da"}}},{"id":"18869148-c21a-990d-387b-4cbc3c4f95ba","type":"GetRegulatoryElementsForGeneInfoFromGene","inputs":{"gene":{"id":"8d325af7-e4a3-00fb-16f8-b26a1c97546d"}}},{"id":"b8a54b2e-77ec-0077-3af2-2045177eef80","type":"GlyGenProteinInformation","inputs":{"gene":{"id":"8d325af7-e4a3-00fb-16f8-b26a1c97546d"}}}],"metadata":{"id":"537e88bb-0fba-fa6c-5586-40c64a6c6545","title":"Use Case 13: Novel Cell Surface Targets for Individual Cancer Patients Analyzed with Common Fund Datasets","description":"The input to this workflow is a data matrix of gene expression that was collected from a pediatric patient tumor patient from the KidsFirst Common Fund program [1]. The RNA-seq samples are the columns of the matrix, and the rows are the raw expression gene count for all human coding genes (Table 1). This data matrix is fed into TargetRanger [2] to screen for targets which are highly expressed in the tumor but lowly expressed across most healthy human tissues based on gene expression data collected from postmortem patients with RNA-seq by the GTEx Common Fund program [3]. Based on this analysis the gene IMP U3 small nucleolar ribonucleoprotein 3 (IMP3) was selected because it was the top candidate returned from the TargetRanger analysis (Tables 2-3). IMP3 is also commonly called insulin-like growth factor 2 mRNA-binding protein 3 (IGF2BP3). Next, we leverage unique knowledge from various other Common Fund programs to examine various functions and knowledge related to IMP3. First, we queried the LINCS L1000 data [4] from the LINCS program [5] converted into RNA-seq-like LINCS L1000 Signatures [6] using the SigCom LINCS API [7] to identify mimicker or reverser small molecules that maximally impact the expression of IMP3 in human cell lines (Fig. 1, Table 4). In addition, we also queried the LINCS L1000 data to identify single gene CRISPR knockouts that down-regulate the expression of IMP3 (Fig. 1, Table 5). These potential drug targets were filtered using the Common Fund IDG program's list of understudied proteins [8] to produce a set of additional targets (Table 6). Next, IMP3 was searched for knowledge provided by the with the Metabolomics Workbench MetGENE tool [9]. MetGENE aggregates knowledge about pathways, reactions, metabolites, and studies from the Metabolomics Workbench Common Fund supported resource [10]. The Metabolomics Workbench was searched to find associated metabolites linked to IMP3 [10]. Furthermore, we leveraged  the Linked Data Hub API [11] to list knowledge about regulatory elements associated with IMP3 (Table 6). Finally, the GlyGen database [12] was queried to identify relevant sets of proteins that are the product of the IMP3 genes, as well as known post-translational modifications discovered on IMP3.\n\nReferences\n1. Heath, Allison P., et al. \"Gabriella Miller Kids First Data Resource Center: Harmonizing clinical and genomic data to support childhood cancer and structural birth defect research.\" Cancer Research 79.13_Supplement (2019): 2464-2464.\n2. G. B. Marino, et al. GeneRanger and TargetRanger: processed gene and protein expression levels across cells and tissues for target discovery, Nucleic Acids Research. 51, W1, 5 July 2023, Pages W213–W224, https://doi.org/10.1093/nar/gkad399\n3. Lonsdale, J. et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics vol. 45 580–585 (2013). doi:10.1038/ng.2653\n4. Subramanian A, et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell. 2017 Nov 30;171(6):1437-1452.e17. doi: 10.1016/j.cell.2017.10.049. \n5. Keenan AB, et al. The Library of Integrated Network-Based Cellular Signatures NIH Program: System-Level Cataloging of Human Cells Response to Perturbations. Cell Syst. 2018 Jan 24;6(1):13-24. doi: 10.1016/j.cels.2017.11.001\n6. Jeon M, et al. Transforming L1000 profiles to RNA-seq-like profiles with deep learning. BMC Bioinformatics. 2022 Sep 13;23(1):374. doi: 10.1186/s12859-022-04895-5.\n7. Evangelista, J. E. et al. SigCom LINCS: data and metadata search engine for a million gene expression signatures. Nucleic Acids Research vol. 50 W697–W709 (2022). doi:10.1093/nar/gkac328\n8. Oprea TI, et al. Unexplored therapeutic opportunities in the human genome. Nat Rev Drug Discov. 2018 May;17(5):317-332. doi: 10.1038/nrd.2018.14\n9. Srinivasan S, et al. MetGENE: gene-centric metabolomics information retrieval tool. Gigascience. 2022 Dec 28;12:giad089. doi: 10.1093/gigascience/giad089\n10. Sud M, et al. Metabolomics Workbench: An international repository for metabolomics data and metadata, metabolite standards, protocols, tutorials and training, and analysis tools. Nucleic Acids Res. 2016 Jan 4;44(D1):D463-70. doi: 10.1093/nar/gkv1042.\n11. Linked Data Hub, https://ldh.genome.network/cfde/ldh/\n12. York, W. S. et al. GlyGen: Computational and Informatics Resources for Glycoscience. Glycobiology vol. 30 72–73 (2019). doi:10.1093/glycob/cwz080\n","summary":"manual","gpt_summary":""}},
  },
]

export default playbooks
