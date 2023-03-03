import * as dict from '@/utils/dict'

export const Disease_backgrounds = [
  {
    name: 'Disease_Signatures_from_GEO_down_2014',
    label: 'Disease Signatures from GEO down',
    termType: 'Disease',
    termRe: /^(?<term>.+) (?<gse>[^ ]+)$/,
    termLabel: 'GEO Disease Signatures',
    termAssociation: 'Containing',
  },
  {
    name: 'Disease_Signatures_from_GEO_up_2014',
    label: 'Disease Signatures from GEO up',
    termType: 'Disease',
    termRe: /^(?<term>.+) (?<gse>[^ ]+)$/,
    termLabel: 'GEO Disease Signatures',
    termAssociation: 'Containing',
  },
  {
    name: 'GTEx_Aging_Signatures_2021',
    label: 'GTEx Aging Signatures 2021',
    termType: 'Disease',
    termRe: /^GTEx (?<term>.+)$/,
    termLabel: 'GeneRIF Rare Diseases',
    termAssociation: 'Containing',
  },
  {
    name: 'Rare_Diseases_GeneRIF_ARCHS4_Predictions',
    label: 'Rare Diseases GeneRIF ARCHS4 Predictions',
    termType: 'Disease',
    termRe: /^(?<term>.+)$/,
    termLabel: 'GeneRIF Rare Diseases',
    termAssociation: 'Containing',
  },
  {
    name: 'Rare_Diseases_GeneRIF_Gene_Lists',
    label: 'Rare Diseases GeneRIF Gene Lists',
    termType: 'Disease',
    termRe: /^(?<term>.+)$/,
    termLabel: 'GeneRIF Rare Diseases',
    termAssociation: 'Containing',
  },
]

export const Drug_backgrounds = [
  {
    name: 'LINCS_L1000_Chem_Pert_Consensus_Sigs',
    label: 'LINCS L1000 Chem Pert Consensus Sigs',
    termType: 'Drug',
    termRe: /^(?<term>.+) (?<direction>Up|Down)$/,
    termLabel: 'L1000 Chem Pert Signatures',
    termAssociation: 'Containing',
  },
  {
    name: 'LINCS_L1000_Chem_Pert_up',
    label: 'LINCS L1000 Chem Pert Up',
    termType: 'Drug',
    termRe: /^(?<desc>.+?-(?<term>.+)-(?<concentration>.+?))$/,
    termLabel: 'L1000 Chem Pert Signatures',
    termAssociation: 'Containing',
  },
  {
    name: 'LINCS_L1000_Chem_Pert_down',
    label: 'LINCS L1000 Chem Pert Down',
    termType: 'Drug',
    termRe: /^(?<desc>.+?-(?<term>.+)-(?<concentration>.+?))$/,
    termLabel: 'L1000 Chem Pert Signatures',
    termAssociation: 'Containing',
  },
]

export const Pathway_backgrounds = [
  {
    name: 'GO_Biological_Process_2021',
    label: 'GO Biological Process 2021',
    termType: 'BiologicalProcess',
    termRe: /^(?<term>.+)\((?<xref>.+?)\)$/,
    termLabel: 'GO Biological Processes',
    termAssociation: 'Containing',
  },
  {
    name: 'KEGG_2019_Human',
    label: 'KEGG 2019 Human',
    termType: 'Pathway',
    termRe: /^(?<term>.+)$/,
    termLabel: 'KEGG Pathways',
    termAssociation: 'Containing',
    extra: {
      organism: 'human',
    },
  },
  {
    name: 'KEGG_2019_Mouse',
    label: 'KEGG 2019 Mouse',
    termType: 'Pathway',
    termRe: /^(?<term>.+)$/,
    termLabel: 'KEGG Pathways',
    termAssociation: 'Containing',
    extra: {
      organism: 'mouse',
    },
  },
  {
    name: 'KEGG_2021_Human',
    label: 'KEGG 2021 Human',
    termType: 'Pathway',
    termRe: /^(?<term>.+)$/,
    termLabel: 'KEGG Pathways',
    termAssociation: 'Containing',
    extra: {
      organism: 'human',
    },
  },
  {
    name: 'MSigDB_Hallmark_2020',
    label: 'MSigDB Hallmark 2020',
    termType: 'Pathway',
    termRe: /^(?<term>.+)$/,
    termLabel: 'MSigDB Hallmark Gene Sets',
    termAssociation: 'Containing',
  },
  {
    name: 'WikiPathway_2021_Human',
    label: 'WikiPathway 2021 Human',
    termType: 'Pathway',
    termRe: /^(?<term>.+)$/,
    termLabel: 'WikiPathways',
    termAssociation: 'Containing',
    extra: {
      organism: 'human',
    },
  },
  {
    name: 'WikiPathway_2019_Human',
    label: 'WikiPathway 2019 Human',
    termType: 'Pathway',
    termRe: /^(?<term>.+) (?<id>[^ ]+)$/,
    termLabel: 'WikiPathways',
    termAssociation: 'Containing',
    extra: {
      organism: 'human',
    },
  },
  {
    name: 'WikiPathway_2019_Mouse',
    label: 'WikiPathway 2019 Mouse',
    termType: 'Pathway',
    termRe: /^(?<term>.+)$/,
    termLabel: 'WikiPathways',
    termAssociation: 'Containing',
    extra: {
      organism: 'mouse',
    },
  },
]

export const Phenotype_backgrounds = [
  {
    name: 'GWAS_Catalog_2019',
    label: 'GWAS Catalog 2019',
    termType: 'Phenotype',
    termRe: /^(?<term>.+)$/,
    termLabel: 'GWAS Phenotypes',
    termAssociation: 'Associated with',
    extra: {
      organism: 'human',
    },
  },
  {
    name: 'UK_Biobank_GWAS_v1',
    label: 'UK Biobank GWAS v1',
    termType: 'Phenotype',
    termRe: /^(?<term>.+)$/,
    termLabel: 'UK Biobank Phenotypes',
    termAssociation: 'Associated with',
    extra: {
      organism: 'human',
    },
  },
  {
    name: 'ClinVar_2019',
    label: 'ClinVar 2019',
    termType: 'Phenotype',
    termRe: /^(?<term>.+)$/,
    termLabel: 'ClinVar Phenotypes',
    termAssociation: 'Associated with',
    extra: {
      organism: 'human',
    },
  },
  {
    name: 'Human_Phenotype_Ontology',
    label: 'Human Phenotype Ontology',
    termType: 'Phenotype',
    termRe: /^(?<term>.+)$/,
    termLabel: 'Human Phenotypes',
    termAssociation: 'Associated with',
    extra: {
      organism: 'human',
    },
  },
  {
    name: 'MGI_Mammalian_Phenotype_Level_4_2019',
    label: 'MGI Mammalian Phenotype Level 4 2019',
    termType: 'Phenotype',
    termRe: /^(?<term>.+)$/,
    termLabel: 'MGI Mammalian Phenotypes',
    termAssociation: 'Associated with',
    extra: {
      organism: 'mouse',
    },
  },
]

export const Tissue_backgrounds = [
  {
    name: 'GTEx_Tissue_Expression_Up',
    label: 'GTEx Tissue Expression Up',
    termType: 'Tissue',
    termRe: /^(?<sample_id>[^ ]+) (?<term>.+?) (?<gender>[^ ]+) (?<age>[^ ]+ years)$/,
    termLabel: 'GTEx Tissue Signatures',
    termAssociation: 'Containing',
    extra: {
      organism: 'human',
    },
  },
  {
    name: 'GTEx_Tissue_Expression_Down',
    label: 'GTEx Tissue Expression Down',
    termType: 'Tissue',
    termRe: /^(?<sample_id>[^ ]+) (?<term>.+?) (?<gender>[^ ]+) (?<age>[^ ]+ years)$/,
    termLabel: 'GTEx Tissue Signatures',
    termAssociation: 'Containing',
    extra: {
      organism: 'human',
    },
  },
  {
    name: 'ARCHS4_Tissues',
    label: 'ARCHS4 Tissues',
    termType: 'Tissue',
    termRe: /^(?<term>.+)$/,
    termLabel: 'ARCHS4 Tissue Signatures',
    termAssociation: 'Containing',
  },
]

export const TranscriptionFactor_backgrounds = [
  {
    name: 'ChEA_2022',
    label: 'ChEA 2022',
    termType: 'TranscriptionFactor',
    termRe: /^(?<term>[^ ]+?) (?<origin>.+) (?<organism>[^ ]+)$/,
    termLabel: 'ChEA Transcription Factors',
    termAssociation: 'Targeting',
  },
  {
    name: 'ENCODE_TF_ChIP-seq_2015',
    label: 'ENCODE TF ChIP-seq 2015',
    termType: 'TranscriptionFactor',
    termRe: /^(?<term>[^ ]+) (?<origin>.+) (?<organism>[^ ]+)$/,
    termLabel: 'ENCODE Transcription Factors',
    termAssociation: 'Targeting',
  },
  {
    name: 'ARCHS4_TFs_Coexp',
    label: 'ARCHS4 TF Co-Expression',
    termType: 'TranscriptionFactor',
    termRe: /^(?<term>[^ ]+) (?<organism>[^ ]+) (?<origin>.+)$/,
    termLabel: 'ARCHS4 Transcription Factors',
    termAssociation: 'Correlated with',
  },
]

export const backgrounds = dict.init([
  ...Disease_backgrounds,
  ...Drug_backgrounds,
  ...Pathway_backgrounds,
  ...Phenotype_backgrounds,
  ...Tissue_backgrounds,
  ...TranscriptionFactor_backgrounds,
].map((value) => ({ key: value.name, value })))
