import type KRG from "@/core/KRG"
import type { FPL } from "@/core/FPPRG"
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { fpl_expand, Metadata, Author } from "./common"
import puppeteer from 'puppeteer';

async function screenshotOf({ graph_id, node_id }: { graph_id: string, node_id: string }) {
  const browser = await puppeteer.launch()
  const page = await browser.newPage()
  await page.goto(`http://localhost:3000/embed/${graph_id}/node/${node_id}`, { waitUntil: 'networkidle0' })
  const pdf = await page.pdf({ format: 'LETTER' })
  await browser.close()
  return pdf
}

function latexEscape(s: string) {
  return s.replaceAll('\\', '\\\\').replaceAll('_', '\\_').replaceAll('{', '\\{').replaceAll('}', '\\}').replaceAll('$', '\\$')
}

export default async function FPL2TEX(props: { krg: KRG, fpl: FPL, metadata?: Metadata, author?: Author | null }): Promise<Record<string, string | Buffer>> {
  const { fullFPL, processLookup, story } = await fpl_expand(props)
  const title = props.metadata?.title ? latexEscape(props.metadata.title) : 'Playbook'
  const abstract = array.unique(story.ast.flatMap(part => !part.tags.includes('abstract') ? [] :
    part.type === 'text' ? [latexEscape(part.text)]
    : part.type === 'cite' ? [`\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' ? [`Fig. \\ref{fig:${story.figures.get(part.ref)}}`]
    : []
  )).join('')
  const introduction = array.unique(story.ast.flatMap(part => !part.tags.includes('introduction') ? [] :
    part.type === 'text' ? [latexEscape(part.text)]
    : part.type === 'cite' ? [`\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' ? [`Fig. \\ref{fig:${story.figures.get(part.ref)}}`]
    : []
  )).join('')
  const methods = array.unique(story.ast.flatMap(part => !part.tags.includes('methods') ? [] :
    part.type === 'text' ? [latexEscape(part.text)]
    : part.type === 'cite' ? [`\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' ? [`Fig. \\ref{fig:${story.figures.get(part.ref)}}`]
    : []
  )).join('')
  const contributors = array.unique(fullFPL.map(head => props.krg.getProcessNode(head.process.type).meta.author).filter((author): author is string => !!author)).map(contributor => latexEscape(contributor))
  // TODO: bibtex for references
  const references = story.ast.flatMap(part => part.type === 'bibitem' ? [`\\bibitem{${story.bibitems.get(part.ref)}}\n${latexEscape(part.text.slice(part.text.indexOf('.')+2))}`] : []).join('\n\n')
  const figures = await Promise.all(fullFPL.map(async (head) => {
    const [figure] = story.ast.filter(part => part.type === 'figure' && part.tags[0] === head.id)
    if (figure?.type !== 'figure') return
    const figure_num = story.figures.get(figure.ref)
    const legend = story.ast.filter(part => part.tags.includes('legend') && part.tags.includes(head.id)).map(part =>
      part.type === 'text' ? latexEscape(part.text)
      : part.type === 'cite' ? `\\cite{${story.bibitems.get(part.ref)}}`
      : part.type === 'figref' ? `Fig. \\ref{fig:${story.figures.get(part.ref)}}`
      : ''
    ).join('')
    return {
      id: head.id,
      files: {
        [`fig${figure_num}.pdf`]: await screenshotOf({ graph_id: fullFPL[fullFPL.length-1].id, node_id: head.id }),
      },
      tex: `
\\begin{figure}[h]
\\centering
\\includegraphics[width=0.9\\textwidth]{fig${figure_num}.pdf}
\\caption{${legend}}\\label{fig:${figure_num}}
\\end{figure}
`,
    }
  }))
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
    ...dict.init(figures.flatMap((fig) => fig ? dict.items(fig.files) : [])),
    'paper.tex': `
\\documentclass{article}
\\usepackage{authblk}
\\usepackage{graphicx}
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

\\section{Figures}\\label{figures}
${figures.flatMap((fig) => fig ? [fig.tex] : []).join('')}

\\clearpage

\\begin{thebibliography}{9}
${references}
\\end{thebibliography}

\\end{document}
`,
    'presentation.tex': `
\\documentclass[11pt]{beamer}
\\usepackage{booktabs}

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
    : part.type === 'cite' ? [`\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' ? [`Fig. \\ref{fig:${story.figures.get(part.ref)}}`]
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
    : part.type === 'cite' ? [`\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' ? [`Fig. \\ref{fig:${story.figures.get(part.ref)}}`]
    : [])).join('')
  const step_figure = figures.flatMap(fig => fig?.id === head.id ? [fig.tex] : []).join('')
  return [
    (step_methods && `
\\subsection{${step_title}}

\\begin{frame}
	\\frametitle{${step_title}}
	
	${step_methods}
\\end{frame}
`) || '',
    (step_figure && `
\\subsubsection{Results}

\\begin{frame}
	\\frametitle{${step_title}: Results}
	
	${step_figure}
\\end{frame}
`) || '',
  ]
}).join('')}

\\begin{frame}
	\\frametitle{References}
	
	\\begin{thebibliography}{99}
    \\footnotesize

    ${references}
  \\end{thebibliography}
\\end{frame}

\\begin{frame}
	\\frametitle{Acknowledgements}
  \\begin{itemize}
    ${contributors.map(contributor => `\\item{${contributor}}`).join('\n')}
  \\end{itemize}
\\end{frame}

\\end{document}
`
  }
}
