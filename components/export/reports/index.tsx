import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { variable_icon } from '@/icons'
import { z } from 'zod'
import { PMCAccessionSet } from '@/components/core/set'
import { GEOAccessionTerm } from '@/components/core/term'
import { PlotlyPlot } from '@/components/viz/plotly'
import { EnrichrScoredGenes, EnrichrScoredPathways, EnrichrScoredPhenotypes } from '@/components/service/enrichr'
import { PerturbSeqrGeneSignature } from '@/components/service/perturbseqr'
import python from '@/utils/python'
import { AnnData } from '@/components/data/anndata'
import { FileC } from '@/components/core/file'
import { useExRouter } from "@/app/fragments/ex-router"


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
        librarySizes: FileC,
        PCAScatter: FileC,
        volcanoScatter: FileC,
        upEnrichrBars: FileC,
        downEnrichrBars: FileC
    }),
    tables: z.object({
        perturbseqrMimickers: FileC,
        perturbseqrReversers: FileC
    }),
    references: z.string().array(),
    model: z.string()
  }))
  .view(props => {
    const router = useExRouter()
    const fpl_id = React.useMemo(() => {
      const m = /^\/(report|graph)\/(.+?)(\/|$)/.exec(router.asPath)
      if (m !== null) return m[1]
    }, [router.asPath])
    return (
      <div className="flex-grow flex flex-col m-0">
        <h2><b>{props.title}</b></h2>
        <h2>Abstract</h2>
        <p>{props.abstract}</p>
        {fpl_id && <>
          <a href={`/api/v1/tex/${fpl_id}?format=pdf`}><button>Download Full Report</button></a>
          <p><i>Can take over a minute to start downloading, please be patient</i></p>
          </>}
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
    enrichrGOBPUp: EnrichrScoredPathways,
    enrichrKEGGUp: EnrichrScoredPathways,
    enrichrChEAUp: EnrichrScoredGenes,
    enrichrKOMPUp: EnrichrScoredPhenotypes,
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
        plots:{
            library_sizes_plot: props.inputs.librarySizesPlot,
            pca_plot: props.inputs.PCAPlot,
            volcano_plot: props.inputs.volcanoPlot
        },
        enrichr_up:{
            enrichr_gobp_up: props.inputs.enrichrGOBPUp,
            enrichr_kegg_up: props.inputs.enrichrKEGGUp,
            enrichr_chea_up: props.inputs.enrichrChEAUp,
            enrichr_komp_up: props.inputs.enrichrKOMPUp
        },
        enrichr_down:{
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