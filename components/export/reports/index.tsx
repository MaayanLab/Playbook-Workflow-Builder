import { MetaNode } from '@/spec/metanode'
import { variable_icon } from '@/icons'
import { z } from 'zod'
import { PMCAccessionSet } from '@/components/core/set'
import { GEOAccessionTerm } from '@/components/core/term'
import { PlotlyPlot } from '@/components/viz/plotly'
import { EnrichrEnrichmentAnalysis, EnrichrScoredGenes, EnrichrScoredPathways, EnrichrScoredPhenotypes } from '@/components/service/enrichr'
import { PerturbSeqrGeneSignature } from '@/components/service/perturbseqr'
import python from '@/utils/python'
import { AnnData } from '@/components/data/anndata'
import { FileC } from '@/components/core/file'
import { ScoredGenes } from '@/components/core/scored'
import { useEffect, useState } from 'react'

const ReportFile = z.object({
  file: FileC,
  caption: z.string()
})

const resolveDRSURL = (fileURL:string) => {
  return fileURL.replace(/^drs:\/\/([^\/]+)\/(.+)$/,
  `${process.env.NEXT_PUBLIC_URL||''}/ga4gh/drs/v1/objects/$2/access/https/data`)
}

export const GEOReanalysisReport = MetaNode(`GEOReanalysisReport`)
  .meta({
    label: 'GEO Re-Analysis Report',
    description: 'An agentic report reanalyzing a GEO study',
    icon: [variable_icon],
    external: true,
  })
  .codec(z.object({
    geo_accession: z.string(),
    title: z.string(),
    abstract: z.string(),
    introduction: z.object({
        problem:z.string(),
        background:z.string(),
        motivation:z.string()
    }),
    methods: z.string(),
    results: z.string(),
    discussion: z.object({
        enrichr: z.string(),
        perturbseqr: z.string(),
        conclusion: z.string()
    }),
    figures: z.object({
        librarySizes: ReportFile,
        PCAScatter: ReportFile,
        volcanoScatter: ReportFile,
        upEnrichrBars: ReportFile,
        downEnrichrBars: ReportFile
    }),
    tables: z.object({
        perturbseqrGeneMimickers: ReportFile,
        perturbseqrDrugMimickers: ReportFile,
        perturbseqrGeneReversers: ReportFile,
        perturbseqrDrugReversers: ReportFile
    }),
    supplement: z.object({
      enrichrUp: z.string(),
      enrichrDown: z.string(),
      perturbseqrUpGenes: z.string(),
      perturbseqrDownGenes: z.string()
    }),
    references: z.array(z.object({
      id: z.string(),
      type: z.string(),
      title: z.string().optional(),
      authors: z.array(z.string()).optional(),
      year: z.string().optional(),
      journal: z.string().optional(),
      volume: z.string().optional(),
      pages: z.string().optional(),
      doi: z.string().optional(),
      url: z.string().optional()
    })),
    model: z.string()
  }))
  .view(props => {
    return (
      <div className="flex-grow flex flex-col m-0 prose">
        <h2><b>{props.title}</b></h2>
        <h2>Abstract</h2>
        <p>{props.abstract}</p>
      </div>
    )
  })
  .build()

export const GEOReanalysisAgentReport = MetaNode(`GEOReanalysisAgentReport`)
  .meta({
    label: 'GEO Re-Analysis Agent Report',
    description: 'Generate an agentic GEO study re-analysis report.',
    icon: [variable_icon],
  })
  .inputs({
    geoAccession: GEOAccessionTerm,
    pmcSet: PMCAccessionSet,
    labelledSamples: AnnData,
    librarySizesPlot: PlotlyPlot,
    PCAPlot: PlotlyPlot,
    volcanoPlot: PlotlyPlot,
    geneSignature: ScoredGenes,
    enrichrUp: EnrichrEnrichmentAnalysis,
    enrichrGOBPUp: EnrichrScoredPathways,
    enrichrKEGGUp: EnrichrScoredPathways,
    enrichrChEAUp: EnrichrScoredGenes,
    enrichrKOMPUp: EnrichrScoredPhenotypes,
    enrichrDown: EnrichrEnrichmentAnalysis,
    enrichrGOBPDown: EnrichrScoredPathways,
    enrichrKEGGDown: EnrichrScoredPathways,
    enrichrChEADown: EnrichrScoredGenes,
    enrichrKOMPDown: EnrichrScoredPhenotypes,
    perturbseqr: PerturbSeqrGeneSignature
  })
  .output(GEOReanalysisReport)
  .resolve(async (props) => {
    return await python(
      'components.export.reports.construct_georeanalysis_report',
      { kargs:[], kwargs: {
        geo_accession: props.inputs.geoAccession,
        pmc_set: props.inputs.pmcSet.set,
        labelled_samples_anndata: props.inputs.labelledSamples,
        signature: props.inputs.geneSignature,
        plots:{
            library_sizes_plot: props.inputs.librarySizesPlot,
            pca_plot: props.inputs.PCAPlot,
            volcano_plot: props.inputs.volcanoPlot
        },
        enrichr_up:{
            enrichr_up_id: props.inputs.enrichrUp,
            enrichr_gobp_up: props.inputs.enrichrGOBPUp,
            enrichr_kegg_up: props.inputs.enrichrKEGGUp,
            enrichr_chea_up: props.inputs.enrichrChEAUp,
            enrichr_komp_up: props.inputs.enrichrKOMPUp
        },
        enrichr_down:{
            enrichr_down_id: props.inputs.enrichrDown,
            enrichr_gobp_down: props.inputs.enrichrGOBPDown,
            enrichr_kegg_down: props.inputs.enrichrKEGGDown,
            enrichr_chea_down: props.inputs.enrichrChEADown,
            enrichr_komp_down: props.inputs.enrichrKOMPDown
        },
        perturbseqr: props.inputs.perturbseqr
      }},
      message => props.notify({ type: 'info', message }),
    )
  })
  .story(props => ({
    abstract: ``,
    methods: ``,
    legend: ``,
  }))
  .build()
