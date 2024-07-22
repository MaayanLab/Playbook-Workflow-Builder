import * as dict from '@/utils/dict'
import {
  archs4_icon,
  glygen_icon,
  gtex_icon,
  idg_icon,
  lincs_icon,
} from '@/icons'
import Citable from '@/utils/citations'

export const Disease_backgrounds = [
  {
    name: 'Disease_Signatures_from_GEO_down_2014',
    label: 'Disease Signatures from GEO down',
    termType: 'Disease',
    termRe: /^(?<term>.+) (?<gse>[^ ]+)$/,
    termLabel: 'GEO Disease Signatures',
    termAssociation: 'Containing',
    ref: Citable.doi('10.1093/nar/gks1193'),
    icon: [],
    tags: {
      'Data Source': {
        GEO: 1,
      }
    },
  },
  {
    name: 'Disease_Signatures_from_GEO_up_2014',
    label: 'Disease Signatures from GEO up',
    termType: 'Disease',
    termRe: /^(?<term>.+) (?<gse>[^ ]+)$/,
    termLabel: 'GEO Disease Signatures',
    termAssociation: 'Containing',
    ref: Citable.doi('10.1093/nar/gks1193'),
    icon: [],
    tags: {
      'Data Source': {
        GEO: 1,
      }
    },
  },
  {
    name: 'GTEx_Aging_Signatures_2021',
    label: 'GTEx Aging Signatures 2021',
    termType: 'Disease',
    termRe: /^GTEx (?<term>.+)$/,
    termLabel: 'GTEx Aging Signatures',
    termAssociation: 'Containing',
    ref: Citable.doi('10.1038/ng.265'),
    icon: [gtex_icon],
    tags: {
      'Data Source': {
        CFDE: 1,
      }
    },
  },
  {
    name: 'Rare_Diseases_GeneRIF_ARCHS4_Predictions',
    label: 'Rare Diseases GeneRIF ARCHS4 Predictions',
    termType: 'Disease',
    termRe: /^(?<term>.+)$/,
    termLabel: 'GeneRIF Rare Diseases',
    termAssociation: 'Containing',
    ref: Citable.cite('GeneRIF, https://www.ncbi.nlm.nih.gov/gene/about-generif'),
    icon: [],
  },
  {
    name: 'Rare_Diseases_GeneRIF_Gene_Lists',
    label: 'Rare Diseases GeneRIF Gene Lists',
    termType: 'Disease',
    termRe: /^(?<term>.+)$/,
    termLabel: 'GeneRIF Rare Diseases',
    termAssociation: 'Containing',
    ref: Citable.cite('GeneRIF, https://www.ncbi.nlm.nih.gov/gene/about-generif'),
    icon: [],
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
    ref: Citable.doi('10.1093/nar/gkac328'),
    icon: [lincs_icon],
    tags: {
      'Data Source': {
        CFDE: 1,
      }
    },
  },
  {
    name: 'LINCS_L1000_Chem_Pert_up',
    label: 'LINCS L1000 Chem Pert Up',
    termType: 'Drug',
    termRe: /^(?<desc>.+?-(?<term>.+)-(?<concentration>.+?))$/,
    termLabel: 'L1000 Chem Pert Signatures',
    termAssociation: 'Containing',
    ref: Citable.doi('10.1016/j.cell.2017.10.049'),
    icon: [lincs_icon],
    tags: {
      'Data Source': {
        CFDE: 1,
      }
    },
  },
  {
    name: 'LINCS_L1000_Chem_Pert_down',
    label: 'LINCS L1000 Chem Pert Down',
    termType: 'Drug',
    termRe: /^(?<desc>.+?-(?<term>.+)-(?<concentration>.+?))$/,
    termLabel: 'L1000 Chem Pert Signatures',
    termAssociation: 'Containing',
    ref: Citable.doi('10.1016/j.cell.2017.10.049'),
    icon: [lincs_icon],
    tags: {
      'Data Source': {
        CFDE: 1,
      }
    },
  },
  {
    name: 'IDG_Drug_Targets_2022',
    label: 'IDG Drug Targets 2022',
    termType: 'Drug',
    termRe: /^(?<term>.+)$/,
    termLabel: 'IDG Drug Target',
    termAssociation: 'Targeting',
    ref: ``,
    icon: [idg_icon],
    tags: {
      'Data Source': {
        CFDE: 1,
      }
    },
  },
]

export const Pathway_backgrounds = [
  {
    name: 'GO_Biological_Process_2023',
    label: 'GO Biological Process 2023',
    termType: 'BiologicalProcess',
    termRe: /^(?<term>.+)\((?<xref>.+?)\)$/,
    termLabel: 'GO Biological Processes',
    termAssociation: 'Containing',
    ref: Citable.doi('10.1038/75556'),
    icon: [],
  },
  {
    name: 'GO_Biological_Process_2021',
    label: 'GO Biological Process 2021',
    termType: 'BiologicalProcess',
    termRe: /^(?<term>.+)\((?<xref>.+?)\)$/,
    termLabel: 'GO Biological Processes',
    termAssociation: 'Containing',
    ref: Citable.doi('10.1038/75556'),
    icon: [],
    hidden: true,
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
    ref: Citable.doi('10.1002/pro.3715'),
    icon: [],
    hidden: true,
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
    ref: Citable.doi('10.1002/pro.3715'),
    icon: [],
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
    ref: Citable.doi('10.1093/nar/gkac963'),
    icon: [],
  },
  {
    name: 'MSigDB_Hallmark_2020',
    label: 'MSigDB Hallmark 2020',
    termType: 'Pathway',
    termRe: /^(?<term>.+)$/,
    termLabel: 'MSigDB Hallmark Gene Sets',
    termAssociation: 'Containing',
    ref: Citable.doi('10.1073/pnas.0506580102'),
    icon: [],
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
    ref: Citable.doi('10.1093/nar/gkaa1024'),
    icon: [],
  },
  {
    name: 'WikiPathways_2019_Human',
    label: 'WikiPathways 2019 Human',
    termType: 'Pathway',
    termRe: /^(?<term>.+) (?<id>[^ ]+)$/,
    termLabel: 'WikiPathways',
    termAssociation: 'Containing',
    extra: {
      organism: 'human',
    },
    ref: Citable.doi('10.1093/nar/gkaa1024'),
    icon: [],
    hidden: true,
  },
  {
    name: 'WikiPathways_2019_Mouse',
    label: 'WikiPathways 2019 Mouse',
    termType: 'Pathway',
    termRe: /^(?<term>.+)$/,
    termLabel: 'WikiPathways',
    termAssociation: 'Containing',
    extra: {
      organism: 'mouse',
    },
    ref: Citable.doi('10.1093/nar/gkaa1024'),
    icon: [],
    hidden: true,
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
    ref: Citable.doi('10.1093/nar/gkac1010'),
    icon: [],
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
    ref: Citable.doi('10.1371/journal.pmed.1001779'),
    icon: [],
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
    ref: Citable.doi('10.1093/nar/gkx1153'),
    icon: [],
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
    ref: Citable.doi('10.1093/nar/gkaa1043'),
    icon: [],
  },
  {
    name: 'MGI_Mammalian_Phenotype_Level_4_2021',
    label: 'MGI Mammalian Phenotype Level 4 2021',
    termType: 'Phenotype',
    termRe: /^(?<term>.+)$/,
    termLabel: 'MGI Mammalian Phenotypes',
    termAssociation: 'Associated with',
    extra: {
      organism: 'mouse',
    },
    ref: Citable.doi('10.1093/nar/gkaa1083'),
    icon: [],
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
    ref: Citable.doi('10.1093/nar/gkaa1083'),
    icon: [],
    hidden: true,
  },
  {
    name: 'KOMP2_Mouse_Phenotypes_2022',
    label: 'KOMP2 Mouse Phenotypes 2022',
    termType: 'Phenotype',
    termRe: /^(?<term>.+?)( \(?<iri>.+?\))?$/,
    termLabel: 'Mouse Phenotypes',
    termAssociation: 'Associated with',
    extra: {
      organism: 'mouse',
    },
    ref: Citable.doi('10.1093/nar/gkac972'),
    icon: [],
    tags: {
      'Data Source': {
        CFDE: 1,
      }
    },
  },
]

export const Tissue_backgrounds = [
  {
    name: 'GTEx_Tissues_V8_2023',
    label: 'GTEx Tissues V8 2023',
    termType: 'Tissue',
    termRe: /^(?<term>.+) (?<gender>[^ ]+) (?<age>[^ ]+) (?<direction>[^ ]+)$/,
    termLabel: 'GTEx Tissue Signatures',
    termAssociation: 'Containing',
    extra: {
      organism: 'human',
    },
    ref: Citable.doi('10.1038/ng.2653'),
    icon: [gtex_icon],
    tags: {
      'Data Source': {
        CFDE: 1,
      }
    },
  },
  {
    name: 'ARCHS4_Tissues',
    label: 'ARCHS4 Tissues',
    termType: 'Tissue',
    termRe: /^(?<term>.+)$/,
    termLabel: 'ARCHS4 Tissue Signatures',
    termAssociation: 'Containing',
    ref: Citable.doi('10.1038/s41467-018-03751-6'),
    icon: [archs4_icon],
    tags: {
      'Data Source': {
        GEO: 1,
      }
    },
  },
  {
    name: 'HuBMAP_ASCTplusB_augmented_2022',
    label: 'HuBMAP ASCTplusB augmented 2022',
    termType: 'Tissue',
    termRe: /^(?<cell_type>.+) - (?<term>.+?)$/,
    termLabel: 'HuBMAP ASCT+B Cell Type Biomarkers',
    termAssociation: 'Containing',
    ref: Citable.cite('HuBMAP ASCT+B Reporter, https://hubmapconsortium.github.io/ccf-asct-reporter/'),
    icon: [],
    tags: {
      'Data Source': {
        CFDE: 1,
      }
    },
  },
  {
    name: 'MoTrPAC_2023',
    label: 'MoTrPAC 2023',
    termType: 'Tissue',
    termRe: /^(\w+?)-(?<term>.+?) (Consensus|(?<gender>Male|Female) (?<timepoint>\w+) (?<dir>Up|Down))$/,
    termLabel: 'MoTrPAC Exercise Signatures',
    termAssociation: 'Containing',
    ref: Citable.doi('10.1016/j.cell.2020.06.004'),
    icon: [],
    tags: {
      'Data Source': {
        CFDE: 1,
      }
    },
  },
]

export const Gene_backgrounds = [
  {
    name: 'ChEA_2022',
    label: 'ChEA 2022',
    termType: 'TranscriptionFactor',
    termRe: /^(?<term>[^ ]+?) (?<origin>.+) (?<organism>[^ ]+)$/,
    termLabel: 'ChEA Transcription Factors',
    termAssociation: 'Targeting',
    ref: Citable.doi('10.1093/nar/gkz446'),
    icon: [],
  },
  {
    name: 'ENCODE_TF_ChIP-seq_2015',
    label: 'ENCODE TF ChIP-seq 2015',
    termType: 'TranscriptionFactor',
    termRe: /^(?<term>[^ ]+) (?<origin>.+) (?<organism>[^ ]+)$/,
    termLabel: 'ENCODE Transcription Factors',
    termAssociation: 'Targeting',
    ref: Citable.doi('10.1038/nature11247'),
    icon: [],
  },
  {
    name: 'ARCHS4_TFs_Coexp',
    label: 'ARCHS4 TF Co-Expression',
    termType: 'TranscriptionFactor',
    termRe: /^(?<term>[^ ]+) (?<organism>[^ ]+) (?<origin>.+)$/,
    termLabel: 'ARCHS4 Transcription Factors',
    termAssociation: 'Correlated with',
    ref: Citable.doi('10.1038/s41467-018-03751-6'),
    icon: [archs4_icon],
  },
  {
    name: 'LINCS_L1000_CRISPR_KO_Consensus_Sigs',
    label: 'LINCS L1000 CRISPR KO Consensus Sigs',
    termType: 'Gene',
    termRe: /^(?<term>.+) (?<direction>.+)$/,
    termLabel: 'L1000 CRISPR KO Signatures',
    termAssociation: 'Containing',
    ref: Citable.doi('10.1093/nar/gkac328'),
    icon: [lincs_icon],
    tags: {
      'Data Source': {
        CFDE: 1,
      }
    },
  },
]

export const Glycan_backgrounds = [
  {
    name: 'GlyGen_Glycosylated_Proteins_2022',
    label: 'GlyGen Glycosylated Proteins 2022',
    termType: 'Glycan',
    termRe: /^(?<term>.+)$/,
    termLabel: 'Glycans',
    termAssociation: 'Glycosylating',
    ref: Citable.doi('10.1093/glycob/cwz080'),
    icon: [glygen_icon],
    tags: {
      'Data Source': {
        CFDE: 1,
      }
    },
  },
]

export const Metabolite_backgrounds = [
  {
    name: 'Metabolomics_Workbench_Metabolites_2022',
    label: 'Metabolomics Workbench Metabolites 2022',
    termType: 'Metabolite',
    termRe: /^(\w+?)-(?<tissue>.+?) (?<gender>Male|Female) (?<timepoint>\w+) (?<dir>Up|Down)$/,
    termLabel: 'Metabolites',
    termAssociation: 'Associated with',
    ref: Citable.cite('The Metabolomics Workbench, https://www.metabolomicsworkbench.org/'),
    icon: [],
    tags: {
      'Data Source': {
        CFDE: 1,
      }
    },
  }
]

export const backgrounds = dict.init([
  ...Disease_backgrounds,
  ...Drug_backgrounds,
  ...Pathway_backgrounds,
  ...Phenotype_backgrounds,
  ...Tissue_backgrounds,
  ...Gene_backgrounds,
  ...Glycan_backgrounds,
  ...Metabolite_backgrounds,
].map((value) => ({ key: value.name, value })))
