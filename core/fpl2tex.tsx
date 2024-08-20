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
  const body = await page.$('body')
  const boundingBox = await body?.boundingBox()
  const pdf = await page.pdf(boundingBox ? { width: boundingBox.width, height: boundingBox.height } : { format: 'letter' })
  await browser.close()
  return pdf
}

function latexEscape(s: string) {
  return s.replaceAll('\\', '\\\\').replaceAll('_', '\\_').replaceAll('{', '\\{').replaceAll('}', '\\}').replaceAll('$', '\\$')
}

export default async function FPL2TEX(props: { krg: KRG, fpl: FPL, metadata?: Metadata, author?: Author | null }): Promise<Record<string, string | Buffer>> {
  const { fullFPL, processLookup, story } = await fpl_expand(props)
  const title = props.metadata?.title ? latexEscape(props.metadata.title) : 'Playbook'
  const abstract = story.ast.flatMap(part => !part.tags.includes('abstract') ? [] :
    part.type === 'text' ? [latexEscape(part.text)]
    : part.type === 'cite' ? [`\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' ? [`Fig. \\ref{fig:${story.figures.get(part.ref)}}`]
    : []
  ).join('')
  const introduction = array.unique(story.ast.flatMap(part => !part.tags.includes('introduction') ? [] :
    part.type === 'text' ? [latexEscape(part.text)]
    : part.type === 'cite' ? [`\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' ? [`Fig. \\ref{fig:${story.figures.get(part.ref)}}`]
    : []
  )).join('')
  const methods = story.ast.flatMap(part => !part.tags.includes('methods') ? [] :
    part.type === 'text' ? [latexEscape(part.text)]
    : part.type === 'cite' ? [`\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' ? [`Fig. \\ref{fig:${story.figures.get(part.ref)}}`]
    : []
  ).join('')
  const contributors = array.unique(fullFPL.map(head => props.krg.getProcessNode(head.process.type).meta.author).filter((author): author is string => !!author)).map(contributor => latexEscape(contributor))
  const figures = (await Promise.all(fullFPL.map(async (head) => {
    const [figure] = story.ast.filter(part => part.type === 'figure' && part.tags[0] === head.id)
    if (figure?.type !== 'figure') return []
    const figure_num = story.figures.get(figure.ref)
    const legend = story.ast.filter(part => part.tags.includes('legend') && part.tags.includes(head.id)).map(part =>
      part.type === 'text' ? latexEscape(part.text)
      : part.type === 'cite' ? `\\cite{${story.bibitems.get(part.ref)}}`
      : part.type === 'figref' ? `Fig. \\ref{fig:${story.figures.get(part.ref)}}`
      : ''
    ).join('')
    return [{
      id: head.id,
      files: {
        [`figures/fig${figure_num}.pdf`]: await screenshotOf({ graph_id: fullFPL[fullFPL.length-1].id, node_id: head.id }),
      },
      label: `fig:${figure_num}`,
      filename: `figures/fig${figure_num}.pdf`,
      legend,
    }]
  }))).flatMap(fig => fig)
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
    ...dict.init(figures.flatMap((fig) => dict.items(fig.files))),
    'paper.tex': `
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

\\section{Figures}\\label{figures}
${figures.flatMap((fig) => fig ? [
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
`,
    'presentation.tex': `
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
	

  ${step_figure.map(fig => `
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
`,
    'poster.tex': `
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

\\documentclass[a0,portrait]{a0poster}

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

${figures.map(fig => `
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
`,
    'references.bib': `
${story.ast.flatMap(part => part.type === 'bibitem' ? part.bibtex ? [`@${part.bibtex.type}{${story.bibitems.get(part.ref)},${part.bibtex.record}}`] : [`@misc{${story.bibitems.get(part.ref)},title={${latexEscape(part.text.slice(part.text.indexOf('.')+2))}}}`] : []).join('\n\n')}
`,
  }
}
