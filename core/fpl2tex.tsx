import type KRG from "@/core/KRG"
import type { FPL } from "@/core/FPPRG"
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { fpl_expand, Metadata, Author } from "./common"
import ReactDOMServer from 'react-dom/server'
import { v4 as uuid4 } from 'uuid'
import fs from 'fs'
import puppeteer from 'puppeteer';

async function screenshotOf(node: React.ReactNode) {
  const html = ReactDOMServer.renderToString(<>{node}</>)
  const tmpdir = `/tmp/${uuid4()}`
  await new Promise<void>((resolve, reject) => fs.mkdir(tmpdir, (err) => {if (err) { reject(err) } else { resolve() }}))
  await fs.writeFileSync(`${tmpdir}/index.html`, html)
  const browser = await puppeteer.launch()
  const page = await browser.newPage()
  await page.goto(`file://${tmpdir}/index.html`, { waitUntil: 'networkidle0' })
  const pdf = await page.pdf({ format: 'LETTER' })
  await browser.close()
  await new Promise<void>((resolve, reject) => fs.unlink(`${tmpdir}/index.html`, (err) => {if (err) { reject(err) } else { resolve() }}))
  await new Promise<void>((resolve, reject) => fs.rmdir(`${tmpdir}`, (err) => {if (err) { reject(err) } else { resolve() }}))
  return pdf
}

export default async function FPL2TEX(props: { krg: KRG, fpl: FPL, metadata?: Metadata, author?: Author | null }): Promise<Record<string, string | Buffer>> {
  const { fullFPL, processLookup, story } = await fpl_expand(props)
  const abstract = story.ast.flatMap(part => !part.tags.includes('abstract') ? [] :
    part.type === 'text' ? [part.text]
    : part.type === 'cite' ? [`\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' ? [`\\ref{fig:${story.figures.get(part.ref)}}`]
    : []
  ).join('')
  const introduction = story.ast.flatMap(part => !part.tags.includes('introduction') ? [] :
    part.type === 'text' ? [part.text]
    : part.type === 'cite' ? [`\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' ? [`\\ref{fig:${story.figures.get(part.ref)}}`]
    : []
  ).join('')
  const methods = story.ast.flatMap(part => !part.tags.includes('methods') ? [] :
    part.type === 'text' ? [part.text]
    : part.type === 'cite' ? [`\\cite{${story.bibitems.get(part.ref)}}`]
    : part.type === 'figref' ? [`\\ref{fig:${story.figures.get(part.ref)}}`]
    : []
  ).join('')
  // TODO: bibtex for references
  const references = story.ast.flatMap(part => part.type === 'bibitem' ? [`\\bibitem{${story.bibitems.get(part.ref)}}\n${part.text.slice(part.text.indexOf('.')+2)}`] : []).join('\n\n')
  const figures = await Promise.all(fullFPL.map(async (head) => {
    const { metanode, output } = processLookup[head.process.id]
    const [figure] = story.ast.filter(part => part.type === 'figure' && part.tags[0] === head.id)
    if (figure?.type !== 'figure') return
    const figure_num = story.figures.get(figure.ref)
    const legend = story.ast.filter(part => part.tags.includes('legend') && part.tags.includes(head.id)).map(part =>
      part.type === 'text' ? part.text
      : part.type === 'cite' ? `\\cite{${story.bibitems.get(part.ref)}}`
      : part.type === 'figref' ? `\\ref{fig:${story.figures.get(part.ref)}}`
      : ''
    ).join('')
    return {
      files: {
        [`fig${figure_num}.pdf`]: await screenshotOf(metanode.output.view(output)),
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
  return {
    ...dict.init(figures.flatMap((fig) => fig ? dict.items(fig.files) : [])),
    'index.tex': `
\\documentclass{article}
\\providecommand{\\keywords}[1]
{
  \\small	
  \\textbf{\\textit{Keywords---}} #1
}

\\begin{document}

${props.metadata?.title ? `\\title{${props.metadata.title}}` : '\\title{Playbook}'}

${props.author ? `\\author${props.author.affiliation ? `[1]` : ''}{${props.author.name}}` : ''}${props.author?.email ? `\\email{${props.author.email}}` : ''}
${props.author?.affiliation ? `\\affil*[1]{${props.author.affiliation}}` : ''}

\\abstract{${abstract}
\\keywords{${[
  'Playbook Workflow Builder',
  ...array.unique(
    dict.values(processLookup)
      .flatMap(({ metanode }) =>
        metanode.meta.tags ? dict.items(metanode.meta.tags).flatMap(({ key: _, value }) => dict.keys(value)).join(' ') : []
      )
  )
].join(', ')}}

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
`
  }
}
