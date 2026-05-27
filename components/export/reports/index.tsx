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

const ReportFigure = z.object({
  pdf: FileC,
  png: FileC,
  caption: z.string()
})

const ReportTable = z.object({
  file: FileC,
  caption: z.string()
})

function renderLatexText(
    text: string, 
    refs: Record<string, {type: string, title: string, id: string, url?:string, authors?: string[], year?: string, journal?: string,volume?: string,pages?: string,doi?: string}>,
    refIndex: Record<string, number>,
    figIndex: Record<string,number>,
    tableIndex: Record<string,number>
  ) {
  return text
    .replace(/~?\\cite\{([^}]+)\}/g, (_, keys) =>
      keys.split(',').map((k: string) => {
        const n = refIndex[k.trim()]
        return n ? `<a href="#ref-${k.trim()}" title="${refs[k.trim()]?.title ?? k}">[${n}]</a>` : `[${k}]`
      }).join('')
    )
    .replace(/\\ref\{fig:([^}]+)\}/g, (_, key) =>
      `<a href="#fig-${key}">${figIndex[key]}</a>`
    )
    .replace(/\\ref\{table:([^}]+)\}/g, (_, key) =>
      `<a href="#fig-${key}">${tableIndex[key]}</a>`
    )
    .replace(/\\href\{([^}]+)\}\{([^}]+)\}/g, (_, link, key) =>
      `<a href="${link}" target="_blank" class="text-blue-600 underline">${key}</a>`
    )
    .replace(/~\\\\/g, '')
}

const Section = ({ title, children }: { title: string, children: React.ReactNode }) => (
  <section className="mt-6">
    <h2 className="text-lg font-bold border-b border-gray-300 mb-2">{title}</h2>
    {children}
  </section>
)

const resolveDRSURL = (fileURL:string) => {
  return fileURL.replace(/^drs:\/\/([^\/]+)\/(.+)$/,
  `${process.env.NEXT_PUBLIC_URL||''}/ga4gh/drs/v1/objects/$2/access/https/data#zoom=page-fit&toolbar=0&navpanes=0&scrollbar=0&zoom=page-width`)
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
        librarySizes: ReportFigure,
        PCAScatter: ReportFigure,
        volcanoScatter: ReportFigure,
        upEnrichrBars: ReportFigure,
        downEnrichrBars: ReportFigure
    }),
    tables: z.object({
        perturbseqrGeneMimickers: ReportTable,
        perturbseqrDrugMimickers: ReportTable,
        perturbseqrGeneReversers: ReportTable,
        perturbseqrDrugReversers: ReportTable
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
      title: z.string(),
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
    const refs = Object.fromEntries(props.references.map(ref => [ref.id, ref]))
    const orderedKeys: string[] = []
    const seenKeys = new Set<string>()
    
    const allTexts = [
      props.abstract,
      props.introduction.problem,
      props.introduction.background,
      props.introduction.motivation,
      props.methods,
      props.results,
      props.discussion.enrichr,
      props.discussion.perturbseqr,
      props.discussion.conclusion,
    ]
    
    const CITE_RE = /\\cite\{([^}]+)\}/g

    for (const text of allTexts) {
      let match: RegExpExecArray | null

      while ((match = CITE_RE.exec(text)) !== null) {
        for (const key of match[1].split(',').map((k: string) => k.trim())) {
          if (key && !seenKeys.has(key)) {
            seenKeys.add(key)
            orderedKeys.push(key)
          }
        }
      }

      CITE_RE.lastIndex = 0
    }
    
    const refIndex = Object.fromEntries(orderedKeys.map((k, i) => [k, i + 1]))

    const figKeyToNum = Object.fromEntries(Object.keys(props.figures).map((k, i) => [k, i + 1]))
    const tableKeyToNum = Object.fromEntries(Object.keys(props.tables).map((k, i) => [k, i + 1]))

    const figureAspectRatios: Record<string, string> = {
      librarySizes: '2 / 1',
      PCAScatter: '13 / 10',
      volcanoScatter: '1 / 1',
      upEnrichrBars: '9 / 5',
      downEnrichrBars: '9 / 5',
    }
  
    const Prose = ({ text }: { text: string }) => (
      <p
        className="my-2 leading-relaxed"
        dangerouslySetInnerHTML={{ __html: renderLatexText(text, refs, refIndex, figKeyToNum, tableKeyToNum) }}
      />
    )

    const [tableData, setTableData] = useState<Record<string, string[][]>>({})
    useEffect(() => {
      async function loadTables() {
        console.log(figKeyToNum)
        const entries = await Promise.all(
          Object.entries(props.tables).map(async ([key, table]) => {
            try {
              const response = await fetch(resolveDRSURL(table.file.url))
              const text = await response.text()

              const parsed = text
                .trim()
                .split("\n")
                .map(line => line.split("\t"))

              return [key, parsed] as const
            } catch (e) {
              console.error(`Failed to load table ${key}`, e)
              return [key, []] as const
            }
          })
        )

        setTableData(Object.fromEntries(entries))
      }

      loadTables()
    }, [props.tables])
  
    return (
      <div className="max-w-6xl mx-auto font-serif px-8 py-6 text-sm">
        <h1 className="text-2xl font-bold text-center mb-4">{props.title}</h1>
  
        <Section title="Abstract"><Prose text={props.abstract} /></Section>
        <Section title="Introduction">
          <Prose text={props.introduction.problem} />
          <Prose text={props.introduction.background} />
          <Prose text={props.introduction.motivation} />
        </Section>
        <Section title="Methods"><Prose text={props.methods} /></Section>
  
        <Section title="Results">
          <Prose text={props.results} />
          <div className="grid grid-cols-1 gap-4 mt-4">
          {Object.entries(props.figures).map(([key, fig]) => (
              <figure key={key} id={`fig-${key}`} className="text-center">
                <img 
                  src={resolveDRSURL(fig.png.url)}
                  style={{
                    width: '75%',
                    //aspectRatio: figureAspectRatios[key],
                    border: 'none',
                    margin: 'auto',
                    display: 'block',
                  }}
                  title={fig.png.filename}
                />
                <figcaption className="text-s mt-1">
                  <b>Figure {figKeyToNum[key]}</b>: <Prose text={fig.caption} />
                </figcaption>
              </figure>
            ))}
          </div>
          <div className="grid grid-cols-1 gap-4 mt-4">
            {Object.entries(props.tables).map(([key, table]) => {
              const rows = tableData[key] || []

              if (rows.length === 0) {
                return (
                  <figure
                    key={key}
                    id={`table-${key}`}
                    className="border p-2"
                  >
                    Loading table...
                  </figure>
                )
              }

              const [header, ...body] = rows

              return (
                <figure
                  key={key}
                  id={`table-${key}`}
                  className="border p-2 overflow-auto"
                >
                  <table className="min-w-full border-collapse text-xs leading-tight">
                    <thead>
                      <tr>
                        {header.map((cell, i) => (
                          <th
                            key={i}
                            className="border px-2 py-1 bg-gray-100 text-left align-top whitespace-nowrap"
                          >
                            {cell}
                          </th>
                        ))}
                      </tr>
                    </thead>

                    <tbody>
                      {body.map((row, i) => (
                        <tr key={i}>
                          {row.map((cell, j) => (
                            <td
                              key={j}
                              className="border px-2 py-1 align-top break-words"
                            >
                              {cell}
                            </td>
                          ))}
                        </tr>
                      ))}
                    </tbody>
                  </table>

                  <figcaption className="text-s mt-1">
                    <b>Table {tableKeyToNum[key]}</b>:{" "}
                    <Prose text={table.caption} />
                  </figcaption>
                </figure>
              )
            })}
          </div>
        </Section>
  
        <Section title="Discussion">
          <Prose text={props.discussion.enrichr} />
          <Prose text={props.discussion.perturbseqr} />
          <Prose text={props.discussion.conclusion} />
        </Section>
      
        <Section title="References">
          <ol className="text-xs space-y-1">
            {orderedKeys.map((key, i) => {
              const ref = refs[key]
              if (!ref) return null

              const authorText =
                ref.authors && ref.authors.length > 0
                  ? ref.authors.length === 1
                    ? ref.authors[0]
                    : `${ref.authors[0]} et al`
                  : ''

              const isWebsite = ref.type==='online'

              return (
                <li key={key} id={`ref-${key}`}>
                  <span className="font-bold">[{i + 1}]</span>{' '}
                  {authorText && <>{authorText}. </>}
                  {isWebsite ? (
                    <>
                      {ref.title}{' '}
                      <a href={ref.url} className="text-blue-600 underline">{ref.url}</a>
                    </>
                  ) : (
                    <>
                      {ref.title}{'. '}
                      {ref.journal && <i>{ref.journal}</i>}
                      {ref.volume && <b>{ref.volume}</b>}
                      {ref.pages && <>{ref.pages}</>}
                      {ref.year && <> ({ref.year})</>}
                      {ref.url && (
                        <>{'. '}<a href={`${ref.url}`} className="text-blue-600 underline">{ref.url}</a></>
                      )}
                    </>
                  )}
                </li>
              )
            })}
          </ol>
        </Section>
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