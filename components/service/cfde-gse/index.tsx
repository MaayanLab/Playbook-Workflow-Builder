import React from "react";
import { MetaNode } from "@/spec/metanode";
import { drc_icon } from "@/icons"
import { z } from 'zod'
import pluralize from "pluralize"
import * as str from '@/utils/str'
import { GeneSet } from "@/components/core/set"
import { Disease, Drug, Gene, Glycan, Metabolite, Phenotype, Tissue } from "@/components/core/primitives"
import { Disease_backgrounds, Drug_backgrounds, Phenotype_backgrounds, Tissue_backgrounds, Gene_backgrounds, Glycan_backgrounds, Metabolite_backgrounds } from "@/components/service/enrichr/backgrounds"
import { EnrichrEnrichmentAnalysis, EnrichrScoredDiseases, EnrichrScoredDrugs, EnrichrScoredGenes, EnrichrScoredGlycans, EnrichrScoredMetabolites, EnrichrScoredPathways, EnrichrScoredPhenotypes, EnrichrScoredTissues, resolveEnrichrGenesetSearchResults } from "@/components/service/enrichr"
import { GMT } from "@/components/data/gene_matrix_transpose";
import python from '@/utils/python'
import SafeRender from "@/utils/saferender";

const cfde_gse_url = 'https://gse.cfde.cloud'
const cfde_gse_enrich_url = 'https://maayanlab.cloud/Enrichr'
const cfde_gse_icon = drc_icon

export const CFDEGSEKG = MetaNode('CFDEGSEKG')
  .meta({
    label: 'CFDE Gene Set Enrichment Knowledge Graph',
    description: 'A CFDE-Centric knowledge graph for gene set enrichment results',
    icon: [cfde_gse_icon],
    external: true,
  })
  .codec(z.union([
    z.object({
      empty: z.literal(true),
    }),
    z.object({
      shortId: z.string(),
      userListId: z.number(),
    })
  ]))
  .view(userlist => {
    const searchParams = new URLSearchParams()
    if (!('empty' in userlist)) {
      searchParams.append('q', JSON.stringify({
        "min_lib": 1,
        "libraries":[
          {"name":"LINCS_L1000_Chem_Pert_Consensus_Sigs","limit":5},
          {"name":"HuBMAP_ASCTplusB_augmented_2022","limit":5},
        ],
        "userListId":userlist.userListId.toString(),
        "search":true,
      }))
      searchParams.append('fullscreen', 'true')
    }
    return (
      <div className="flex-grow flex flex-row m-0" style={{ minHeight: 1000 }}>
        {'empty' in userlist ? 
          <div className="prose max-w-none">Enrichment Analysis Cannot be Performed on Empty Set</div>
          : <iframe
            className="flex-grow border-0"
            src={`${cfde_gse_url}?${searchParams.toString()}`}
          />}
      </div>
    )
  })
  .build()

export const EnrichrUserListToCFDEGSEKG = MetaNode('EnrichrUserListToCFDEGSEKG')
  .meta({
    label: 'Visualize CFDE Gene Set Enrichment Knowledge Graph',
    description: 'View gene set results as a knowledge graph',
    tags: {
      'Output Type': {
        'Knowledge Graph': 1,
      },
      'Data Source': {
        'CFDE': 1
      },
    },
    external: true,
  })
  .inputs({ enrichr: EnrichrEnrichmentAnalysis })
  .output(CFDEGSEKG)
  .resolve(async (props) => props.inputs.enrichr)
  .story(props => ({
    introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
    methods: `The enrichment results from Enrichr\\ref{doi:10.1002/cpz1.90} are used to subset a Knowledge Graph of inter-connected CFDE Gene Set Libraries using CFDE GSE\\ref{CFDE GSE, https://gse.cfde.cloud/}.`,
    legend: `An interactive page containing a knowledge graph highlighting the enriched terms in the CFDE Gene Set Libraries provided by CFDE GSE\\ref{CFDE GSE, https://gse.cfde.cloud/}.`,
  }))
  .build()

export const CFDEGSEGenesetSearch = MetaNode('CFDEGSEGenesetSearch')
  .meta({
    label: `CFDE Enrichment Analysis`,
    tags: {
      'Input Type': {
        Gene: 1,
      },
      'Input Cardinality': {
        Set: 1,
      },
      Service: {
        'CFDE GSE': 1,
      },
      'Output Type': {
        Signature: 1,
        Gene: 1,
      },
      'Output Cardinality': {
        Weighted: 1,
      },
    },
    icon: [cfde_gse_icon],
    description: 'Perform Enrichment Analysis',
    external: true,
  })
  .inputs({ geneset: GeneSet })
  .output(CFDEGSEKG)
  .resolve(async (props) => {
    if (props.inputs.geneset.set.length === 0) {
      return { empty: true }
    }
    const formData = new FormData()
    formData.append('list', props.inputs.geneset.set.join('\n'))
    formData.append('description', str.truncate(`playbook-partnership${props.inputs.geneset.description ? `:${props.inputs.geneset.description}` : ''}`, 255))
    const response = await fetch(`${cfde_gse_enrich_url}/addList`, {
      method: 'post',
      body: formData,
    })
    if (response.status !== 200) {
      throw new Error(`'CFDE GSE addList' Request failed with status ${response.status}: ${await response.text()}`)
    }
    return await response.json()
  })
  .story(props => ({
    abstract: `The gene set${props.inputs && props.inputs.geneset.description ? ` containing ${props.inputs.geneset.description}` : ''} was submitted to CFDE GSE\\ref{Common Fund Data Ecosystem (CFDE) Gene Set Enrichment (GSE), https://gse.cfde.cloud}.`,
    introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. CFDE GSE is a gene set search engine & knowledge graph that enables the querying of CF annotated gene sets\\ref{Common Fund Data Ecosystem (CFDE) Gene Set Enrichment (GSE), https://gse.cfde.cloud}.`,
    methods: `The gene set in ${props.input_refs?.geneset} is submitted to CFDE GSE\\ref{Common Fund Data Ecosystem (CFDE) Gene Set Enrichment (GSE), https://gse.cfde.cloud} for enrichment analysis.`,
    legend: `An interactive page with the results of the CFDE GSE\\ref{Common Fund Data Ecosystem (CFDE) Gene Set Enrichment (GSE), https://gse.cfde.cloud} enrichment analysis results.`,
  }))
  .build()

export const CFDEGSEGenesetSearchT = [
  { backgrounds: Disease_backgrounds.filter(bg => 'tags' in bg && bg.tags && 'Data Source' in bg.tags && bg.tags['Data Source']['CFDE'] === 1),
    output: EnrichrScoredDiseases, T: Disease },
  { backgrounds: Drug_backgrounds.filter(bg => 'tags' in bg && bg.tags && 'Data Source' in bg.tags && bg.tags['Data Source']['CFDE'] === 1),
    output: EnrichrScoredDrugs, T: Drug },
  { backgrounds: Phenotype_backgrounds.filter(bg => 'tags' in bg && bg.tags && 'Data Source' in bg.tags && bg.tags['Data Source']['CFDE'] === 1),
    output: EnrichrScoredPhenotypes, T: Phenotype },
  { backgrounds: Tissue_backgrounds.filter(bg => 'tags' in bg && bg.tags && 'Data Source' in bg.tags && bg.tags['Data Source']['CFDE'] === 1),
    output: EnrichrScoredTissues, T: Tissue },
  { backgrounds: Gene_backgrounds.filter(bg => 'tags' in bg && bg.tags && 'Data Source' in bg.tags && bg.tags['Data Source']['CFDE'] === 1),
    output: EnrichrScoredGenes, T: Gene },
  { backgrounds: Glycan_backgrounds.filter(bg => 'tags' in bg && bg.tags && 'Data Source' in bg.tags && bg.tags['Data Source']['CFDE'] === 1),
    output: EnrichrScoredGlycans, T: Glycan },
  { backgrounds: Metabolite_backgrounds.filter(bg => 'tags' in bg && bg.tags && 'Data Source' in bg.tags && bg.tags['Data Source']['CFDE'] === 1),
    output: EnrichrScoredMetabolites, T: Metabolite },
].flatMap(({ backgrounds, output, T }) =>
backgrounds.map(bg => ({ bg, output, T }))
).map(({ bg, output, T }) =>
  MetaNode(`ExtractCFDEGSEGenesetSearch[${bg.name}]`)
    .meta({
      label: `Extract Enriched ${bg.termLabel}`,
      icon: [cfde_gse_icon, ...(bg.icon||[])],
      description: `Extract Significant Terms from the ${bg.label} Library`,
      tags: {
        'Output Type': {
          [pluralize(T.label)]: 1,
        },
        ...('tags' in bg ? bg.tags as Record<string, Record<string, 1>> : {}),
      },
      external: true,
    })
    .inputs({ searchResults: CFDEGSEKG })
    .output(output)
    .resolve(async (props) => {
      return !('empty' in props.inputs.searchResults) ? {
        background: bg.name,
        ...(await resolveEnrichrGenesetSearchResults(bg, props.inputs.searchResults))
      } : {
        background: bg.name,
        terms: [],
        scored: [],
      }
    })
    .story(props => ({
      abstract: `The gene set was enriched against the ${bg.label}${bg.ref} library to identify statistically significant ${bg.termLabel}.`,
      introduction: `Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets\\ref{doi:10.1002/cpz1.90}.`,
      methods: `The enrichment analysis results against the ${bg.label}${bg.ref} library are resolved via the Enrichr API\\ref{doi:10.1002/cpz1.90}.`,
      tableLegend: `A table of significantly enriched ${bg.termLabel} from the ${bg.label}${bg.ref} library paired with their significance scores as reported by Enrichr\\ref{doi:10.1002/cpz1.90}.`,
    }))
    .build()
)

const datasets = {
  "GlyGen":{
      "name":"GlyGen Glycosylated Proteins",
      "resource":"GlyGen",
      "url":"https://cfde-drc.s3.amazonaws.com/GlyGen/XMT/2022-12-13/GlyGen_XMT_2022-12-13_GlyGen_Glycosylated_Proteins_2022.gmt",
      "processing":"The GlyGen Glycosylated Proteins dataset was processed from human proteoform glycosylation site citation datasets from GlyGen ~\\cite{GlyGen}. Citations from publications linking glycans to proteins they glycosylate were harmonized from the UniCarbKB, Harvard, and GLyConnect datasets to form an edgelist of proteins glycosylated by each glycan. Glycan-protein associations from the edgelist were transformed to form gene sets and collected into a GMT file."
  },
  "GTExAging":{
      "name":"GTEx Tissue-Specific Aging Signatures",
      "resource":"GTEx",
      "url":"https://cfde-drc.s3.amazonaws.com/GTEx/XMT/2022-06-06/GTEx_XMT_2022-06-06_GTEx_Aging_Signatures_2021.gmt",
      "processing":"The GTEx Tissue-Specific Aging dataset was constructed using gene expression count data from the Genotype Tissue Expression (GTEx) Portal ~\\cite{GTEx}. Sample metadata was filtered to create groups of samples with matching tissue type and donor age ranges. Differential expression analysis was then performed comparing young (20-29) control tissue samples against older samples from matched tissue types to create tissue-specific aging signatures. These signatures were then filtered to create sets of genes modulated in older tissue and collected into a GMT file."
  },
  "GTExTissues":{
    "name":"GTEx Tissues",
    "resource":"GTEx",
    "url":"https://cfde-drc.s3.amazonaws.com/GTEx/XMT/2023-03-10/GTEx_XMT_2023-03-10_GTEx_Tissues_V8_2023.gmt",
    "processing":"The GTEx Tissues dataset was constructed using gene expression count data from  the Genotype Tissue Expression (GTEx) Portal. Sample metadata was filtered to create groups of samples with matching tissue type and donor age ranges. The grouped data was log2-transformed, quantile normalized, and z-scored to create up and down gene sets."
  },
  "HuBMAPASCTB":{
    "name":"HuBMAP ASCT+B",
    "resource":"HuBMAP",
    "url":"https://cfde-drc.s3.amazonaws.com/HuBMAP/XMT/2025-03-21/HuBMAP_ASCT-B.gmt",
    "processing":"The HuBMAP ASCT+B Augmented dataset expands sets of biomarker genes for anatomical structures and cell types. A GMT of known biomarkers was constructed using the ASCT+B table for each tissue from HuBMAP. Tissue-specific gene correlation matrices were constructed by filtering for RNA-seq expression from ARCHS4 from each tissue, performing TMM and log normalization, and calculating pairwise gene correlation. These tissue-specific matrices were used to retrieve the 100 most genes most correlated to the original ASCT+B genes, creating augmented AS/CT gene sets."
  },
  "IDGDrugTargets":{
    "name":"IDG Drug Targets",
    "resource":"IDG",
    "url":"https://cfde-drc.s3.amazonaws.com/IDG/XMT/2022-12-13/IDG_XMT_2022-12-13_IDG_Drug_Targets_2022.gmt",
    "processing":"The IDG Drug Targets dataset was constructed using drug-target interaction data from DrugCentral. For each drug, the associated targets were collected into a gene set, keeping gene sets with a length of at least 5."
  },
  "KOMP2":{
    "name":"KOMP2 Mouse Phenotypes",
    "resource":"KOMP2",
    "url":"https://cfde-drc.s3.amazonaws.com/KOMP2/XMT/2022-12-13/KOMP2_XMT_2022-12-13_KOMP2_Mouse_Phenotypes_2022.gmt",
    "processing":"The KOMP2 Mouse Phenotypes dataset was constructed using knockout mouse genotype-phenotype data from KOMP2. Gene sets were converted to human orthologs and sets with at least 5 genes were kept."
  },
  "LINCSChemical":{
      "name":"LINCS L1000 Chemical Perturbation Signatures",
      "resource":"LINCS",
      "url":"https://s3.perturbseqr.maayanlab.cloud/data/lincs-l1000-cp.gmt.gz",
      "processing":"The LINCS L1000 Chemical Perturbation Signatures dataset was processed from Level 3 L1000 data from the CLUE platform ~\\cite{LINCS}. Replicate perturbation profiles were grouped by matching timepoint, perturbagen, dosage, detection plate, and well ID metadata. Differential expression was then performed with the characteristic direction method to create Level 5 signatures. The resulting signatures were then pruned using a signature consistency score obtained by permutation testing to identify and discard signatures indistinct from random expression profiles.  Remaining signatures were then collected into a GMT file."
  },
  "LINCSCRISPR":{
      "name":"LINCS L1000 CRISPR KO Perturbation Signatures",
      "resource":"LINCS",
      "url":"https://s3.perturbseqr.maayanlab.cloud/data/lincs-l1000-xpr.gmt.gz",
      "processing":"The LINCS L1000 CRISPR KO Perturbation Signatures dataset was processed from Level 3 L1000 data from the CLUE platform ~\\cite{LINCS}. Replicate perturbation profiles were grouped by matching timepoint, perturbagen, detection plate, and well ID metadata. Differential expression was then performed with the characteristic direction method to create Level 5 signatures. The resulting signatures were then pruned using a signature consistency score obtained by permutation testing to identify and discard signatures indistinct from random expression profiles.  Remaining signatures were then collected into a GMT file."
  },
  "Metabolomics":{
      "name": "Metabolomics Metabolite Enzyme Associations",
      "resource":"Metabolomics",
      "url":"https://cfde-drc.s3.amazonaws.com/MW/XMT/2022-12-13/MW_XMT_2022-12-13_Metabolomics_Workbench_Metabolites_2022.gmt",
      "processing":"The Metabolomics Enzyme-Metabolite Associations dataset was processed from associations between metabolites and related enzymes ~\\cite{Metabolomics}. The Metabolomics Workbench metabolite dataset was created by integrating structures and annotations from various open resources. The Human Metabolome Gene/Protein (MGP) database containing information about metabolites and linked genes/proteins was used to create an edgelist of metabolite-enzyme associations for metabolites with at least one public study deposited to the Metabolomics Workbench. The ensyme-metabolite edgelist was harmonized to create gene sets which were collected into a GMT."
  },
  "MoTrPAC":{
      "name":"MoTrPAC Rat Endurance Exercise Training Transcriptomics",
      "resource":"MoTrPAC",
      "url":"https://cfde-drc.s3.amazonaws.com/MoTrPAC/XMT/2024-03-05/MoTrPAC_Endurance_Trained_Rats_2023.gmt",
      "processing":"The MoTrPAC Rat Endurance Exercise Training dataset was processed from transcriptomic data from the Molecular Transducers of Physical Activity Consortium (MoTrPAC) Data Hub ~\\cite{MoTrPAC}. Differential expression data from MoTrPAC was filtered to keep associations where gene expression changed over the course of 8 weeks of endurance exercise. Exercise tissue samples were differentiated by sex and tissue type in addition to exercise timepoint. The downloaded data also included training data of genes with modulated expression over the training window in addition to the timewise data. The differential expression data was filtered to create gene sets and a GMT."
  },
  "RummaGEODrug":{
      "name":"RummaGEO Drug Perturbation Signatures",
      "resource":"RummaGEO",
      "url":"https://cfde-drc.s3.amazonaws.com/LINCS/XMT/2025-06-25/RummaGEO_DrugPerturbations_2025.gmt",
      "processing":"RummaGEO is a webserver developed to enable expression signature search against human and mouse RNA-seq studies deposited to GEO ~\\cite{RummaGEO}. To enable this search, offline processing automatically identified groups of samples were  based on shared conditions and differential expression analysis was performed. The RummaGEO Drug Perturbation Signatures dataset was processed by querying the RummaGEO metadata search API to identify a subset of signatures where one condition included a control term and the other included a drug. Gene sets for each signature in the subset were combined to construct a GMT file."
  },
  "RummaGEOGene":{
      "name":"RummaGEO Gene Perturbation Signatures",
      "resource":"RummaGEO",
      "url":"https://cfde-drc.s3.amazonaws.com/LINCS/XMT/2025-06-25/RummaGEO_GenePerturbationSignatures.gmt.txt",
      "processing":"RummaGEO is a webserver developed to enable expression signature search against human and mouse RNA-seq studies deposited to GEO ~\\cite{RummaGEO}. To enable this search, offline processing automatically identified groups of samples were  based on shared conditions and differential expression analysis was performed. The RummaGEO Gene Perturbation Signatures dataset was processed by querying the RummaGEO metadata search API to identify a subset of signatures where one condition included a control term and the other included a gene perturbation. Gene sets for each signature in the subset were combined to construct a GMT file."
  }
}

export const CFDEDataset = MetaNode('CFDEDataset')
.meta({
  label: 'CFDE Dataset',
  description: 'A dataset from the CFDE Workbench',
    tags: {
      'Data Source': {
        'CFDE': 1
      },
    },
    external: true,
  })
  .codec(z.object({
    key:z.string(),
    dataset: z.object({
      name: z.string(),
      resource: z.string(),
      url: z.string(),
      processing: z.string()
    })
  }))
  .view(dataset => {
    return (<></>)
  })
  .build()

export const CFDEGMTs = MetaNode('CFDEGMTs')
.meta({
  label: 'CFDE GMTs',
  description: 'Multiple GMTs from the CFDE Workbench',
    tags: {
      'Data Source': {
        'CFDE': 1
      },
    },
    external: true,
  })
  .codec(z.array(z.object({
    key:z.string(),
    dataset: z.object({
      name: z.string(),
      resource: z.string(),
      url: z.string(),
      processing: z.string(),
      gmt: z.record(
        z.string(), 
        z.object({ 
          description: z.string().optional(), 
          set: z.array(z.string()) 
        })
      )
    })
  })))
  .view(datasets => {
    return (
      <></>
    )
  })
  .build()

export const SelectCFDEGSEDataset = MetaNode('SelectCFDEGSEDataset')
.meta({
  label: 'Select a CFDE GSE Dataset',
  description: 'Select a dataset from the CFDE Workbench',
    tags: {
      'Data Source': {
        'CFDE': 1
      },
    },
    external: true,
  })
  .codec(z.object({ key: z.string() }))
  .inputs({})
  .output(CFDEDataset)
  .prompt(props => {
    const [selected, setSelected] = React.useState<string | null>(
      props.output?.key ?? null
    )

    return (
      <>
        <div className="md-6">
          <select
            className="select select-bordered w-full"
            value={selected ?? ''}
            onChange={e => setSelected(e.target.value)}
          >
            <option value="" disabled>Select a dataset...</option>
            {Object.entries(datasets).map(([key, dataset]) => (
              <option key={key} value={key}>{dataset.name}</option>
            ))}
          </select>
        </div>
        <button
          className="bp5-button bp5-large"
          disabled={selected === null}
          onClick={() => {
            if (selected !== null) props.submit({ key: selected })
          }}
        >
          Submit
        </button>
        {props.output ? <SafeRender component={CFDEDataset.view} props={props.output} /> : null}
      </>
    )
  })
  .resolve(async (props) => {
    return {key: props.data.key, dataset: datasets[props.data.key as keyof typeof datasets]}
  })
  .story(props => ({
    abstract: `${props.output ? 'The ' + props.output.dataset.name : 'A'} dataset was selected.`,
    methods: `${props.output ? 'The ' + props.output.dataset.name : 'A'} dataset was selected.`,
    }))
  .build()

export const CFDEGSEDatasetToGMT = MetaNode('CFDEGSEDatasetToGMT')
  .meta({
    label: 'CFDE GSE Dataset To GMT',
    description: 'Retrieve a GMT for a dataset from the CFDE Workbench',
      tags: {
        'Data Source': {
          'CFDE': 1
        },
      },
      external: true,
    })
    .inputs({dataset:CFDEDataset})
    .output(GMT)
    .resolve(async (props) => {
      const {dataset} = props.inputs.dataset
      const res = await fetch(dataset.url)
      if (!res.ok) throw new Error(`Failed to fetch GMT from ${dataset.url}`)
      return await python(
        'components.service.cfde-gse.load_gmt',
        { kargs: [], kwargs: { url:dataset.url } },
        message => props.notify({ type: 'info', message }),
      )
    })
    .story(props => ({
      abstract: `${props.inputs ? 'The ' + props.inputs.dataset.dataset.name : 'A'} GMT was retrieved from the CFDE Workbench\\ref{doi:10.1016/j.jmb.2026.169631}.`,
    }))
    .build()

export const SelectCFDEGSEGMTs = MetaNode('SelectCFDEGSEGMTs')
  .meta({
    label: 'Select CFDE GSE GMTs',
    description: 'Select multiple datasets from the CFDE Workbench and retrieve GMTs',
      tags: {
        'Data Source': {
          'CFDE': 1
        },
      },
      external: true,
    })
    .codec(z.object({ 
      datasets: z.array(z.string()),
    }))
    .inputs({})
    .output(CFDEGMTs)
    .prompt(props => {
      const [selected0, setSelected0] = React.useState<string | null>(
        props.data?.datasets[0] ?? null
      )
      const [selected1, setSelected1] = React.useState<string | null>(
        props.data?.datasets[1] ?? null
      )
      const [selected2, setSelected2] = React.useState<string | null>(
        props.data?.datasets[2] ?? null
      )
      const [selected3, setSelected3] = React.useState<string | null>(
        props.data?.datasets[3] ?? null
      )
      const [selected4, setSelected4] = React.useState<string | null>(
        props.data?.datasets[4] ?? null
      )
  
      return (
        <>
          <div className="md-6">
            <select
              className="select select-bordered w-full"
              value={selected0 ?? ''}
              onChange={e => setSelected0(e.target.value)}
            >
              <option value="" disabled>Select a dataset...</option>
              {Object.entries(datasets).map(([key, dataset]) => (
                <option key={key} value={key}>{dataset.name}</option>
              ))}
            </select>
            <select
              className="select select-bordered w-full"
              value={selected1 ?? ''}
              onChange={e => setSelected1(e.target.value)}
            >
              <option value="" disabled>Select a dataset...</option>
              {Object.entries(datasets).map(([key, dataset]) => (
                <option key={key} value={key}>{dataset.name}</option>
              ))}
            </select>
            <select
              className="select select-bordered w-full"
              value={selected2 ?? ''}
              onChange={e => setSelected2(e.target.value)}
            >
              <option value="" disabled>Select a dataset...</option>
              {Object.entries(datasets).map(([key, dataset]) => (
                <option key={key} value={key}>{dataset.name}</option>
              ))}
            </select>
            {selected0 && selected1 && selected2 && 
              <select
                className="select select-bordered w-full"
                value={selected3 ?? ''}
                onChange={e => setSelected3(e.target.value)}
              >
                <option value="" disabled>Select a dataset...</option>
                {Object.entries(datasets).map(([key, dataset]) => (
                  <option key={key} value={key}>{dataset.name}</option>
                ))}
              </select>
            }
            {selected0 && selected1 && selected2 && selected3 && 
              <select
                className="select select-bordered w-full"
                value={selected4 ?? ''}
                onChange={e => setSelected4(e.target.value)}
              >
                <option value="" disabled>Select a dataset...</option>
                {Object.entries(datasets).map(([key, dataset]) => (
                  <option key={key} value={key}>{dataset.name}</option>
                ))}
              </select>
            }
          </div>
          <button
            className="bp5-button bp5-large"
            disabled={!selected0 || !selected1 || !selected2}
            onClick={() => {
              if (!selected0 || !selected1 || !selected2) {
                throw new Error("Please select at least 3 datasets.")
              }
              props.submit({ datasets:[selected0, selected1, selected2, selected3, selected4].filter((ds): ds is string => ds !== null) })
            }}
          >
            Submit
          </button>
          {props.output ? <SafeRender component={CFDEGMTs.view} props={props.output} /> : null}
        </>
      )
    })
    .resolve(async (props) => {
      console.log(props.data)
      const selected = Object.entries(datasets).filter(([key]) => props.data.datasets.includes(key))
      return await python(
        'components.service.cfde-gse.load_gmts',
        { kargs: [], kwargs: { datasets:selected } },
        message => props.notify({ type: 'info', message }),
      )
    })
    .story(props => ({
      abstract: `${props.output ? `The ${props.output .map(d => d.dataset.name).slice(0, -1).join(", ")}, and ${props.output.at(-1)!.dataset.name} datasets were selected from the CFDE Workbench\\ref{doi:10.1016/j.jmb.2026.169631}.`: ""}`
      }))
    .build()