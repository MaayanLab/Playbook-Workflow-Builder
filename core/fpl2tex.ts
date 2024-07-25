import type KRG from "@/core/KRG"
import type { FPL } from "@/core/FPPRG"
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { fpl_expand, Metadata, Author } from "./common"

function screenshotOf(node: React.ReactNode) {
  return 'TODO'
}

export default async function FPL2TEX(props: { krg: KRG, fpl: FPL, metadata?: Metadata, author?: Author | null }): Promise<string> {
  const { fullFPL, processLookup, story } = await fpl_expand(props)
  const abstract = story.ast.flatMap(part => !part.tags.includes('abstract') ? [] :
    part.type === 'text' ? [part.text]
    : part.type === 'cite' ? [`\\cite{${part.ref}}`]
    : part.type === 'ref' ? [`\\ref{${part.ref}}`]
    : []
  ).join('')
  const introduction = story.ast.flatMap(part => !part.tags.includes('introduction') ? [] :
    part.type === 'text' ? [part.text]
    : part.type === 'cite' ? [`\\cite{${part.ref}}`]
    : part.type === 'ref' ? [`\\ref{${part.ref}}`]
    : []
  ).join('')
  const methods = story.ast.flatMap(part => !part.tags.includes('methods') ? [] :
    part.type === 'text' ? [part.text]
    : part.type === 'cite' ? [`\\cite{${part.ref}}`]
    : part.type === 'ref' ? [`\\ref{${part.ref}}`]
    : []
  ).join('')
  const references = story.ast.flatMap(part => part.type === 'bibitem' ? [`\\bibitem{${part.ref}}\n${part.text.slice(part.text.indexOf('.')+2)}`] : []).join('\n\n')
  return `
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
${
  fullFPL
    .map((head) => {
      const { metanode, output } = processLookup[head.process.id]
      const [figure] = story.ast.filter(part => part.type === 'figure' && part.tags.includes(head.process.id))
      if (figure?.type !== 'figure') return ''
      const legend = story.ast.filter(part => part.tags.includes('legend') && part.tags.includes(head.process.id)).map(part =>
        part.type === 'text' ? part.text
        : part.type === 'cite' ? `\\cite{${part.ref}}`
        : part.type === 'ref' ? `\\ref{${part.ref}}`
        : ''
      ).join('')
      return `
\\begin{figure}[h]
\\centering
\\includegraphics[width=0.9\\textwidth]{${screenshotOf(metanode.output.view(output))}}
\\caption{${legend}}\\label{${figure.ref}}
\\end{figure}
`
  }).join('')
}

\\clearpage

\\begin{thebibliography}{9}
${references}
\\end{thebibliography}

\\end{document}
`
}
