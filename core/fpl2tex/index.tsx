import type KRG from "@/core/KRG"
import type { Data, FPL } from "@/core/FPPRG"
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { fpl_expand, Metadata, Author } from "@/core/common"
import puppeteer from 'puppeteer'
import path from 'path'
import fs from 'fs'
import cache from '@/utils/global_cache'

const puppeteerSingleton = cache('puppeteer', () => puppeteer.launch())

type ReportFile = {
  url: string,
  filename: string,
  description?: string,
  size: number,
  sha256: string
} 

type ReanalysisReport = {
  geo_accession: string,
  title: string,
  abstract: string,
  introduction: {
    problem: string,
    background: string,
    motivation: string
  },
  methods: string,
  results: string,
  discussion: {
    enrichr: string,
    perturbseqr: string,
    conclusion: string
  },
  figures : {
    librarySizes: ReportFile,
    PCAScatter: ReportFile,
    volcanoScatter: ReportFile,
    upEnrichrBars: ReportFile,
    downEnrichrBars: ReportFile
  },
  tables: {
    perturbseqrMimickers: ReportFile,
    perturbseqrReversers: ReportFile
  },
  references: string[]
}

async function screenshotOf({ graph_id, node_id }: { graph_id: string, node_id: string }) {
  const browser = await puppeteerSingleton
  const ctx = await browser.createBrowserContext()
  const page = await ctx.newPage()
  await page.goto(`http://localhost:3000/embed/${graph_id}/node/${node_id}`)
  await page.waitForSelector('div.embed-ready')
  await page.waitForNetworkIdle({ idleTime: 500 })
  const body = await page.$('body')
  const boundingBox = await body?.boundingBox()
  const pdf = await page.pdf(boundingBox ? { width: boundingBox.width, height: boundingBox.height, format: 'letter' } : { width: 2551, height: 3295, format: 'letter' })
  await ctx.close()
  return pdf
}

function latexEscape(s: string) {
  return s.replace(/(\\|_|\{|\}|\$|&)/g, '\\$1')
}

async function extras() {
  const extrasRootPath = path.resolve(
    process.env.APP_ROOT as string,
    'core',
    'fpl2tex',
    'extras',
  )
  return {
    'extras/a0poster.cls': new Promise<ArrayBuffer>((resolve, reject) =>
      fs.readFile(path.resolve(extrasRootPath, 'a0poster.cls'), (err, data) => {
        if (err) { reject(err) } else { resolve(new Uint8Array(data).buffer) }
      })
    ),
    'extras/a0size.sty': new Promise<ArrayBuffer>((resolve, reject) =>
      fs.readFile(path.resolve(extrasRootPath, 'a0size.sty'), (err, data) => {
        if (err) { reject(err) } else { resolve(new Uint8Array(data).buffer) }
      })
    ),
  }
}

export async function fetchGEOReportFile( url:string ): Promise<string | ArrayBuffer> {
  const objecturl = url.replace(
    /^drs:\/\/([^\/]+)\/(.+)$/,
    `${process.env.NEXT_PUBLIC_URL||''}/ga4gh/drs/v1/objects/$2/access/https/data`
  )
  const response = await fetch(objecturl);

  if (!response.ok) {
    throw new Error(`Failed to fetch: ${response.status} ${response.statusText}`);
  }

  return response.arrayBuffer();
}

export async function constructGEOReportReferences( references:string[] ): Promise<string> {
  return references.join('\n\n')
}

export async function GEOReanalysis2TEX(report_id:string,{ geo_accession, title, abstract,introduction,methods,results,discussion,figures,tables,references }:ReanalysisReport): Promise<Record<string, Promise<string | ArrayBuffer>>> {
  async function GEOextras() {
    const extrasRootPath = path.resolve(
      process.env.APP_ROOT as string,
      'core',
      'fpl2tex',
      'extras',
    )
    
    return {
      'ExcelAtFIT.cls': new Promise<ArrayBuffer>((resolve, reject) =>
        fs.readFile(path.resolve(extrasRootPath, 'ExcelAtFIT.cls'), (err, data) => {
          if (err) { reject(err) } else { resolve(new Uint8Array(data).buffer) }
        })
      ),
      'maayan-lab-logo.png': new Promise<ArrayBuffer>((resolve, reject) =>
        fs.readFile(path.resolve(extrasRootPath, 'maayan-lab-logo.png'), (err, data) => {
          if (err) { reject(err) } else { resolve(new Uint8Array(data).buffer) }
        })
      )
    }
  }
  const texFile = `${geo_accession}.tex`

  const files:Record<string, Promise<string | ArrayBuffer>> = {
    ...(await GEOextras()),
    'figures/librarySizeBars.pdf': fetchGEOReportFile(figures.librarySizes.url),
    'figures/PCAScatter.pdf': fetchGEOReportFile(figures.PCAScatter.url),
    'figures/VolcanoScatter.pdf': fetchGEOReportFile(figures.volcanoScatter.url),
    'figures/UpEnrichrBars.pdf': fetchGEOReportFile(figures.upEnrichrBars.url),
    'figures/DownEnrichrBars.pdf': fetchGEOReportFile(figures.downEnrichrBars.url),
    'tables/PerturbseqrMimickers.csv': fetchGEOReportFile(tables.perturbseqrMimickers.url),
    'tables/PerturbseqrReversers.csv': fetchGEOReportFile(tables.perturbseqrReversers.url),
    [ texFile ]: Promise.resolve(`
  \\documentclass{ExcelAtFIT}
  \\ExcelFinalCopy

  \\hypersetup{
    pdftitle={Paper Title},
    pdfauthor={Author},
    pdfkeywords={Keywords}
  }

  \\lstset{ 
    backgroundcolor=\\color{white},
    basicstyle=\\footnotesize\\tt,
  }

  \\usepackage{float}
  \\usepackage{csquotes}
  \\usepackage{pgfplotstable}
  \\pgfplotsset{compat=1.18}
  \\usepackage{adjustbox}

  \\addbibresource{references.bib}
  
  \\providecommand{\\keywords}[1]
  {
    \\small	
    \\textbf{\\textit{Keywords---}} #1
  }
  
  \\ExcelYear{2026}

  \\PaperTitle{${title}}
  
  \\Authors{Axiom K. Playwright (AI Author)*}
  \\affiliation{*%
    {The Ma'ayan Laboratory, Mount Sinai Center for Bioinformatics, Department of Pharmacological Sciences, Windreich Department of Artificial Intelligence and Human Health, Icahn School of Medicine at Mount Sinai, New York, NY 10029}}

  \\Abstract{
    ${abstract}
    The executed playbook is available at \\url{${process.env.NEXT_PUBLIC_URL||''}/report/${report_id}}.
  }
  
  \\begin{document}


  \\startdocument
  \\raggedbottom

  \\section{Introduction}\\label{introduction}
  ${introduction.problem}
  
  ${introduction.background}

  ${introduction.motivation}


  \\section{Methods}\\label{methods}

  The workflow starts with selecting ${geo_accession} as the search term.
  
  GEO studies were identified matching ${geo_accession} using  ARCHS4 ~\\cite{ARCHS4} term search.
  
  The GEO study accession was used to fetch the linked publication accession from PMC.
  
  Gene expression counts and sample metadata for published samples were obtained from ARCHS4 ~\\cite{ARCHS4}.
  
  An AnnData file was prepared from the input data and metadata ~\\cite{AnnData}.
  
  Genes from the anndata matrix were filtered to include protein-coding genes.
  
  The samples were then labeled as either control or perturbation to allow for further analysis.

  \\begin{figure}[!htbp]
    \\centering
    \\includegraphics[width=\\linewidth]{figures/librarySizeBars.pdf}
    \\caption{A bar plot representing library sizes.}
    \\label{fig:Libraries}
  \\end{figure}
  
  The AnnData file was then visualized as a bar plot representing library sizes (Figure \\ref{fig:Libraries}).
  
  Dimensionality reduction of the data was performed using PCA with the normalization set to log-counts-per-million (logCPM).

  \\begin{figure}[htbp]
    \\centering
    \\includegraphics[width=\\linewidth]{figures/PCAScatter.pdf}
    \\caption{A scatterplot showing the first two principal components of the AnnData.}
    \\label{fig:PCA}
  \\end{figure}
  
  The first two principal components (PCs) were used to generate a scatter plot (Figure \\ref{fig:PCA}).
  
  The AnnData file was then analyzed using differential expression by Limma-Voom ~\\cite{limma, voom} to create a gene signature using the selected conditions. 

  \\begin{figure}[!htbp]
    \\centering
    \\includegraphics[width=\\linewidth]{figures/VolcanoScatter.pdf}
    \\caption{A volcano plot showing the log2-fold-changes and statistical significance of each gene.}
    \\label{fig:Volcano}
  \\end{figure}
  
  The data in the differential expression table was then visualized as a volcano plot (Figure \\ref{fig:Volcano}).
  
  The up-regulated genes were extracted from the gene signature computed by the Limma-Voom analysis from the file.
  
  The gene set containing significant up genes was extracted from the gene signature and submitted to Enrichr ~\\cite{Enrichr}.

  \\begin{figure*}[!htbp]
    \\centering
    \\includegraphics[width=\\linewidth]{figures/UpEnrichrBars.pdf}
    \\caption{Results of the Enrichr ~\\cite{Enrichr} enrichment analysis for down-regulated genes. Bar charts show the significantly enriched terms from the GO Biological Process 2023 ~\\cite{GO}, KEGG 2021 Human ~\\cite{KEGG}, ChEA 2022 ~\\cite{ChEA}, and KOMP2 Mouse Phenotypes 2022 ~\\cite{KOMP2} libraries. Z-scores are capped at 10x the smallest value shown.}
    \\label{fig:Enrichr_up}
  \\end{figure*}
  
  The gene set was enriched against the GO Biological Process 2023 ~\\cite{GO}, KEGG 2021 Human ~\\cite{KEGG}, ChEA 2022 ~\\cite{ChEA}, and KOMP2 Mouse Phenotypes 2022 ~\\cite{KOMP2} libraries to identify statistically significant enriched biological processes, pathways, transcription factors and phenotypes (Figure \\ref{fig:Enrichr_up}). 
  
  The gene set containing significant down genes was extracted from the gene signature and submitted to Enrichr ~\\cite{Enrichr}.
  
  \\begin{figure*}[!htbp]
    \\centering
    \\includegraphics[width=\\linewidth]{figures/DownEnrichrBars.pdf}
    \\caption{Results of the Enrichr ~\\cite{Enrichr} enrichment analysis for down-regulated genes. Bar charts show the significantly enriched terms from the GO Biological Process 2023 ~\\cite{GO}, KEGG 2021 Human ~\\cite{KEGG}, ChEA 2022 ~\\cite{ChEA}, and KOMP2 Mouse Phenotypes 2022 ~\\cite{KOMP2} libraries. Z-scores are capped at 10x the smallest value shown.}
    \\label{fig:Enrichr_down}
  \\end{figure*}

  The gene set was enriched against the GO Biological Process 2023 ~\\cite{GO}, KEGG 2021 Human ~\\cite{KEGG}, ChEA 2022 ~\\cite{ChEA}, and KOMP2 Mouse Phenotypes 2022 ~\\cite{KOMP2} libraries to identify statistically significant enriched biological processes, pathways, transcription factors and phenotypes (Figure \\ref{fig:Enrichr_down}).
  
  \\catcode\`\\#=12
  \\pgfplotstableset{
    every column/.style={string type, column type=l},
    every cell/.style={
        postproc cell content/.append code={
          \\pgfkeysalso{@cell content=\\expandafter\\detokenize\\expandafter{\\pgfkeysvalueof{/pgfplots/table/@cell content}}}
        }
    },
    every head row/.style={before row=\\toprule, after row=\\midrule},
    every last row/.style={after row=\\bottomrule}
  }

  \\begin{table*}[!htbp]
    \\centering
    \\caption{A listing of small molecules and single gene CRISPR KOs that produce gene expression profiles similar to the signatures found using Perturb-seqr ~\\cite{Perturb-Seqr}.}
    \\adjustbox{max width=\\textwidth}{
      \\large
      \\pgfplotstabletypeset[col sep=comma]{tables/PerturbseqrMimickers.csv}
    }
    \\label{table:PerturbSeqrMimic}
  \\end{table*}

  \\begin{table*}[!htbp]
    \\centering
    \\caption{A listing of small molecules and single gene CRISPR KOs that produce gene expression profiles opposite to the signatures found using Perturb-seqr ~\\cite{Perturb-Seqr}.}
    \\adjustbox{max width=\\textwidth}{
      \\large
      \\pgfplotstabletypeset[col sep=comma]{tables/PerturbseqrReversers.csv}
    }
    \\label{table:PerturbSeqrReverse}
  \\end{table*}
  \\catcode\`\\#=6

  Significant genes were extracted from the gene signature  and submitted to Perturb-Seqr ~\\cite{Perturb-Seqr} to identify small molecules and single gene CRISPR KOs producing gene expression profiles similar (Table ~\\ref{table:PerturbSeqrMimic}) or opposite (Table ~\\ref{table:PerturbSeqrReverse}) to the signature.
  
  \\section{Results}\\label{results}

  ${results}
  
  \\section{Discussion}\\label{Discussion}
  ${discussion.enrichr}

  ${discussion.perturbseqr}

  ${discussion.conclusion}
  
  \\section*{Acknowledgements}\\label{Acknowledgements}
  This project was funded by the NIH grant \\href{https://reporter.nih.gov/search/SdeFoZSP2U2zRTjMZKFHlQ/project-details/11080094}{OT2OD036435}.
  The workflow was constructed using the \\href{https://playbook-workflow-builder.cloud/}{Playbook Workflow Builder} and full details about the site and its development are available in its PLOS Computational Biology \\href{https://doi.org/10.1371/journal.pcbi.1012901}{publication}.
  The OpenAI gpt-5-mini LLM was used to assist in writing the introduction.
  
  {\\small
    \\printbibliography
  }
  
  \\end{document}
  `),'references.bib': constructGEOReportReferences(references)
  }

  return files
}

//export default async function FPL2TEX(props: { krg: KRG, fpl: FPL, metadata?: Metadata, author?: Author | null }): Promise<Record<string, string | Buffer>>
export default async function FPL2TEX(props: { krg: KRG, fpl: FPL, metadata?: Metadata, author?: Author | null }): Promise<Record<string, Promise<string | ArrayBuffer>>> {
  const { fullFPL, processLookup, story } = await fpl_expand(props)
  // Use GEO report, if exists
  const reportNode = fullFPL.at(-1)
  if (reportNode && reportNode.process?.type === 'GEOReanalysisAgentReport') {
    const reportSections:ReanalysisReport = (await reportNode.process.output() as Data).value
    return GEOReanalysis2TEX(props.fpl.id, reportSections)
  }
  const title = props.metadata?.title ? latexEscape(props.metadata.title) : 'Playbook'
  const abstract = story.ast.flatMap(part => !part.tags.includes('abstract') ? [] :
    part.type === 'text' ? [latexEscape(part.text)]
    : part.type === 'cite' ? [`~\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'figure' ? [`Fig. \\ref{fig:${story.figures.get(part.ref)?.ref}}`]
    : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'table' ? [`Table \\ref{fig:${story.figures.get(part.ref)?.ref}}`]
    : []
  ).join('')
  const introduction = array.unique(story.ast.flatMap(part => !part.tags.includes('introduction') ? [] :
    part.type === 'text' ? [latexEscape(part.text)]
    : part.type === 'cite' ? [`~\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'figure' ? [`Fig. \\ref{fig:${story.figures.get(part.ref)?.ref}}`]
    : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'table' ? [`Table \\ref{fig:${story.figures.get(part.ref)?.ref}}`]
    : []
  )).join('')
  const methods = story.ast.flatMap(part => !part.tags.includes('methods') ? [] :
    part.type === 'text' ? [latexEscape(part.text)]
    : part.type === 'cite' ? [`~\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'figure' ? [`Fig. \\ref{fig:${story.figures.get(part.ref)?.ref}}`]
    : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'table' ? [`Table \\ref{fig:${story.figures.get(part.ref)?.ref}}`]
    : []
  ).join('')
  const contributors = array.unique(fullFPL.map(head => props.krg.getProcessNode(head.process.type).meta.author).filter((author): author is string => !!author)).map(contributor => latexEscape(contributor))
  const figures = (await array.mapSequential(fullFPL, async (head) => {
    const [figure] = story.ast.filter(part => part.type === 'figure' && part.tags[0] === head.id)
    if (figure?.type !== 'figure') return []
    const figure_num = story.figures.get(figure.ref)?.ref
    const legend = story.ast.filter(part => (part.tags.includes('legend') || part.tags.includes('figureLegend') || part.tags.includes('tableLegend')) && part.tags.includes(head.id)).map(part =>
      part.type === 'text' ? latexEscape(part.text)
      : part.type === 'cite' ? `~\\cite{${story.bibitems.get(part.ref)}}`
      : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'figure' ? `Fig. \\ref{fig:${story.figures.get(part.ref)?.ref}}`
      : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'table' ? `Table \\ref{fig:${story.figures.get(part.ref)?.ref}}`
      : ''
    ).join('')
    return [{
      id: head.id,
      kind: figure.kind,
      files: {
        [`${figure.kind}s/${figure_num}.pdf`]: screenshotOf({ graph_id: fullFPL[fullFPL.length-1].id, node_id: head.id }),
      },
      label: `fig:${figure_num}`,
      filename: `${figure.kind}s/${figure_num}.pdf`,
      legend,
    }]
  })).flatMap(fig => fig)
  const keywords = [
    'Playbook Workflow Builder',
    ...array.unique(
      dict.values(processLookup)
        .flatMap(({ metanode }) =>
          metanode.meta.tags ? dict.items(metanode.meta.tags).flatMap(({ key: _, value }) => dict.keys(value)).map(tag => latexEscape(tag)).join(' ') : []
        )
    )
  ].join(', ')
  
  return {
    ...(await extras()),
    ...dict.init(figures.flatMap((fig) => dict.items(fig.files))),
    'paper.tex': Promise.resolve(`
\\documentclass{article}
\\usepackage[T1]{fontenc}
\\usepackage{authblk}
\\usepackage{graphicx}

\\usepackage[backend=biber,style=numeric,citestyle=nature]{biblatex}
\\addbibresource{references.bib}

\\providecommand{\\keywords}[1]
{
  \\small	
  \\textbf{\\textit{Keywords---}} #1
}

\\begin{document}

\\title{${title}}

${props.author ? `\\author${props.author.affiliation ? `[1]` : ''}{${latexEscape(props.author.name)}}` : ''}${props.author?.email ? ` \\\\ ${latexEscape(props.author.email)}` : ''}
${props.author?.affiliation ? `\\affil[1]{${latexEscape(props.author.affiliation)}}` : ''}

\\abstract{${abstract}}
\\keywords{${keywords}}

\\maketitle

\\section{Introduction}\\label{introduction}
${introduction}

\\section{Methods}\\label{methods}
${methods}

\\section{Results}\\label{results}

\\section{Conclusion}\\label{conclusion}

\\clearpage
\\section{Tables}\\label{tables}
${figures.filter(fig => fig.kind === 'table').flatMap((fig) => fig ? [
  `
\\begin{table}[h]
\\centering
\\includegraphics[width=0.9\\textwidth]{${fig.filename}}
\\caption{${fig.legend}}\\label{${fig.label}}
\\end{table}
`
] : []).join('')}

\\clearpage
\\section{Figures}\\label{figures}
${figures.filter(fig => fig.kind === 'figure').flatMap((fig) => fig ? [
  `
\\begin{figure}[h]
\\centering
\\includegraphics[width=0.9\\textwidth]{${fig.filename}}
\\caption{${fig.legend}}\\label{${fig.label}}
\\end{figure}
`
] : []).join('')}

\\clearpage

\\printbibliography

\\end{document}
`),
    'presentation.tex': Promise.resolve(`
\\documentclass[11pt]{beamer}
\\usepackage[T1]{fontenc}
\\usepackage{booktabs}

\\usepackage[backend=biber,style=numeric,citestyle=nature]{biblatex}
\\addbibresource{references.bib}

% see beamer for available themes
\\usetheme{default}
\\usecolortheme{default}
\\usefonttheme{default}
\\useinnertheme{default}
\\useoutertheme{default}

\\title{${title}}

${props.author ? `\\author${props.author.affiliation ? `[1]` : ''}{${latexEscape(props.author.name)}}` : ''}${props.author?.email ? ` \\\\ ${latexEscape(props.author.email)}` : ''}

${props.author?.affiliation ? `\\institute{${latexEscape(props.author.affiliation)}}` : ''}

\\date{\\today}

\\begin{document}

\\begin{frame}
  \\titlepage
\\end{frame}

\\begin{frame}
  \\frametitle{Presentation Overview}
  \\tableofcontents
\\end{frame}

\\section{Introduction}
${fullFPL.flatMap(head => {
  const process = props.krg.getProcessNode(head.process.type)
  const step_title = latexEscape(process.meta.label)
  const step_introduction = array.unique(story.ast.flatMap(part => !part.tags.includes('introduction') || !part.tags.includes(head.id) ? [] :
    part.type === 'text' ? [latexEscape(part.text)]
    : part.type === 'cite' ? [`~\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'figure' ? [`Fig. \\ref{fig:${story.figures.get(part.ref)?.ref}}`]
    : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'table' ? [`Table \\ref{fig:${story.figures.get(part.ref)?.ref}}`]
    : [])).join('')
  if (!step_introduction) return []
  else return [
    `
\\subsection{${step_title}}

\\begin{frame}
	\\frametitle{${step_title}}
	
	${step_introduction}
\\end{frame}
`
  ]
}).join('')}

\\section{Workflow}

${fullFPL.flatMap(head => {
  const process = props.krg.getProcessNode(head.process.type)
  const step_title = latexEscape(process.meta.label)
  const step_methods = array.unique(story.ast.flatMap(part => !part.tags.includes('methods') || !part.tags.includes(head.id) ? [] :
    part.type === 'text' ? [latexEscape(part.text)]
    : part.type === 'cite' ? [`~\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'figure' ? [`Fig. \\ref{fig:${story.figures.get(part.ref)?.ref}}`]
    : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'table' ? [`Table \\ref{fig:${story.figures.get(part.ref)?.ref}}`]
    : [])).join('')
  const step_figure = figures.filter(fig => fig?.id === head.id)
  return [
    (step_methods && `
\\subsection{${step_title}}

\\begin{frame}
	\\frametitle{${step_title}}
	
	${step_methods}
\\end{frame}
`) || '',
  (step_figure.length > 0 && `
\\subsubsection{Results}

\\begin{frame}
	\\frametitle{${step_title}: Results}
	
  ${step_figure.filter(fig => fig.kind === 'table').map(fig => `
  \\begin{table}[h]
  \\centering
  \\includegraphics[width=0.5\\textwidth]{${fig.filename}}
  \\caption{${fig.legend}}\\label{${fig.label}}
  \\end{table}
`).join('')}

  ${step_figure.filter(fig => fig.kind === 'figure').map(fig => `
  \\begin{figure}[h]
  \\centering
  \\includegraphics[width=0.5\\textwidth]{${fig.filename}}
  \\caption{${fig.legend}}\\label{${fig.label}}
  \\end{figure}
`).join('')}
\\end{frame}
`) || '',
]
}).join('')}

\\begin{frame}
	\\frametitle{References}
	\\printbibliography[heading=none]
\\end{frame}

\\begin{frame}
	\\frametitle{Acknowledgements}
  \\begin{itemize}
    ${contributors.map(contributor => `\\item{${contributor}}`).join('\n')}
  \\end{itemize}
\\end{frame}

\\end{document}
`),
    'poster.tex': Promise.resolve(`
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a0poster Portrait Poster
% LaTeX Template
% Version 1.0 (22/06/13)
%
% The a0poster class was created by:
% Gerlinde Kettl and Matthias Weiser (tex@kettl.de)
% 
% This template has been downloaded from:
% http://www.LaTeXTemplates.com
%
% License:
% CC BY-NC-SA 3.0 (http://creativecommons.org/licenses/by-nc-sa/3.0/)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\\documentclass[a0,portrait]{extras/a0poster}

\\usepackage{multicol} % This is so we can have multiple columns of text side-by-side
\\columnsep=100pt % This is the amount of white space between the columns in the poster
\\columnseprule=3pt % This is the thickness of the black line between the columns in the poster

\\usepackage[svgnames]{xcolor} % Specify colors by their 'svgnames', for a full list of all colors available see here: http://www.latextemplates.com/svgnames-colors

\\usepackage{times} % Use the times font
%\\usepackage{palatino} % Uncomment to use the Palatino font

\\usepackage{graphicx} % Required for including images
\\usepackage{booktabs} % Top and bottom rules for table
\\usepackage[font=small,labelfont=bf]{caption} % Required for specifying captions to tables and figures
\\usepackage{amsfonts, amsmath, amsthm, amssymb} % For math fonts, symbols and environments
\\usepackage{wrapfig} % Allows wrapping text around tables and figures

\\begin{document}

\\begin{minipage}[b]{0.75\\linewidth}
\\veryHuge \\color{NavyBlue} \\textbf{${title}} \\color{Black}\\\\ % Title
\\Huge\\textit{A Playbook Workflow}\\\\[2cm] % Subtitle
${props.author ? `\\huge \\textbf{${latexEscape(props.author.name)}}\\\\[0.5cm] % Author(s)` : ''}
${props.author?.affiliation ? `\\huge {${latexEscape(props.author.affiliation)}}\\\\[0.4cm] % University/organization` : ''}
${props.author?.email ? `\\Large \\texttt{${latexEscape(props.author.email)}}\\\\` : ''}
\\end{minipage}

%\\begin{minipage}[b]{0.25\\linewidth}
%\\includegraphics[width=20cm]{logo.png}\\\\
%\\end{minipage}

\\vspace{1cm} % A bit of extra whitespace between the header and poster content

\\begin{multicols}{2} % This is how many columns your poster will be broken into, a portrait poster is generally split into 2 columns

\\color{Navy} % Navy color for the abstract

\\begin{abstract}

${abstract}

\\end{abstract}

\\color{SaddleBrown} % SaddleBrown color for the introduction

\\section*{Introduction}

${introduction}

\\color{DarkSlateGray} % DarkSlateGray color for the rest of the content

\\section*{Methods}

${methods}

\\section*{Results}

${figures.filter(fig => fig.kind === 'table').map(table => `
\\begin{center}\\vspace{1cm}
\\includegraphics[width=0.5\\linewidth]{${table.filename}}
\\captionof{table}{${table.legend}}\\label{${table.label}}
\\end{center}\\vspace{1cm}
`).join('\n\n')}

${figures.filter(fig => fig.kind === 'figure').map(fig => `
\\begin{center}\\vspace{1cm}
\\includegraphics[width=0.5\\linewidth]{${fig.filename}}
\\captionof{figure}{${fig.legend}}\\label{${fig.label}}
\\end{center}\\vspace{1cm}
`).join('\n\n')}

\\color{SaddleBrown} % SaddleBrown color for the conclusions to make them stand out

\\section*{Conclusions}

% TODO

\\color{DarkSlateGray} % Set the color back to DarkSlateGray for the rest of the content

\\section*{Future Directions}

% TODO

\\bibliographystyle{naturemag}
\\bibliography{references}

\\section*{Acknowledgements}

% TODO

\\end{multicols}
\\end{document}
`),
    'references.bib': Promise.resolve(`
${story.ast.flatMap(part => part.type === 'bibitem' ? part.bibtex ? [`@${part.bibtex.type}{${story.bibitems.get(part.ref)},${part.bibtex.record}}`] : [`@misc{${story.bibitems.get(part.ref)},title={${latexEscape(part.text.slice(part.text.indexOf('.')+2))}}}`] : []).join('\n\n')}
`),
  }
}
