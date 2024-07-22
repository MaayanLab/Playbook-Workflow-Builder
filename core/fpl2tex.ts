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
  return `
\\documentclass{article}
\\begin{document}

${props.metadata?.title ? `\\title{${props.metadata.title}}` : ''}

${props.author ? `\\author${props.author.affiliation ? `[1]` : ''}{${props.author.name}}` : ''}${props.author?.email ? `\\email{${props.author.email}}` : ''}
${props.author?.affiliation ? `\\affil*[1]{${props.author.affiliation}}` : ''}

\\abstract{${story['abstract'].flatMap(token =>
  token.type === 'text' ? [token.text]
  : token.type === 'cite' ? [`\\cite${token.ref}`]
  : []
).join(' ')}}
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

\\section{Introduction}\\label{sec1}
${story['introduction'].flatMap(token =>
  token.type === 'text' ? [token.text]
  : token.type === 'cite' ? [`\\cite${token.ref}`]
  : []
).join('')}

\\section{Methods}\\label{sec2}
${story['methods'].flatMap(token =>
  token.type === 'text' ? [token.text]
  : token.type === 'cite' ? [`\\cite${token.ref}`]
  : []
).join('')}

\\section{Results}\\label{sec3}

\\section{Conclusion}\\label{sec4}

\\section{Figures}\\label{sec5}
${
  dict.values(processLookup)
    .map(({ metanode, output, story }, i) => `
\\begin{figure}[h]
\\centering
\\includegraphics[width=0.9\\textwidth]{${screenshotOf(metanode.output.view(output))}}
\\caption{${story?.legend}}\\label{fig${i+1}}
\\end{figure}
`).join('')
}

\\begin{thebibliography}
${story['references'].flatMap(token =>
  token.type === 'ref' ? [`\\bibitem{${token.ref}} ${token.text}`]
    : []
).join('\n')}
\\end{thebibliography}

\\end{document}
`
}
