import * as array from '@/utils/array'
import * as dict from '@/utils/dict'
import doiMappings from './citations'

function convertDOI(text: string) {
  const m = /^doi:(.+)$/.exec(text)
  if (m !== null && m[1] in doiMappings) {
    const citation = doiMappings[m[1] as keyof typeof doiMappings]
    return { text: citation.nature, bibtex: citation.bibtex, doi: text }
  } else if (m !== null) {
    console.warn(`${m[1]} not registered, call 'npm run codegen:citations'`)
  }
  return { text }
}

/**
 * This lets you add citations inline and pulls them out and puts them at the end.
 *
 * Usage:
 * ```
 * Using some thing\ref{Some Citation},\ref{Some Other Citation}... And also\ref{Some Citation}.
 * ```
 * Output:
 * ```
 * Using some thing [1,2]... And also [1].
 *
 * 1. Some Citation
 * 2. Some Other Citation
 * ```
 * 
 * TODO: make the logic behind the way figrefs are getting determined easier to follow
 *  part of the issue is that figrefs can appear in different order as the figures
 *  so we can't provide a figref number till after the figures are inserted.
 */
export default function extractCitations(texts: { text?: string, tags: string[] }[]) {
  const ast: (
    { type: 'text', text: string, tags: string[] }
    | { type: 'cite', text: string, ref: string, tags: string[] }
    | { type: 'bibitem', text: string, ref: string, bibtex?: { type: string, record: string }, tags: string[] }
    | { type: 'figref', text: string, ref: string, tags: string[] }
    | { type: 'figure', text: string, ref: string, tags: string[] }
  )[] = []
  const figures = new Map<string, string>()
  const figure_asts: Record<string, { type: 'figure', ref: string, text: string, tags: string[] }> = {}
  const bibitems = new Map<string, string>()
  const bibitem_asts: Record<string, { type: 'bibitem', ref: string, text: string, bibtex?: { type: string, record: string }, tags: string[] }> = {}
  for (const { text, tags } of texts) {
    if (!text) continue
    if (tags.includes('legend')) {
      const ref = tags.filter(tag => tag !== 'legend').join('-')
      if (!(ref in figure_asts)) {
        figure_asts[ref] = { type: 'figure', text: '', ref, tags }
      } else {
        figure_asts[ref].tags = array.unique([...figure_asts[ref].tags, ...tags])
      }
      if (!figures.has(ref)) {
        figure_asts[ref].text = ''
        figures.set(ref, `${figures.size+1}`)
        ast.push(figure_asts[ref])
      }
    }
    const expr = /\\(figref|ref)\{(.+?)\}/g
    let i = 0
    let m
    while ((m = expr.exec(text)) !== null) {
      const reftype = m[1]
      const ref = m[2]
      if (reftype === 'ref') {
        if (!bibitems.has(ref)) {
          const bibitem = convertDOI(ref)
          bibitem_asts[ref] = { type: 'bibitem', text: [`${bibitems.size+1}.`, bibitem.text, bibitem.doi].filter((text): text is string => !!text).join(' '), bibtex: bibitem.bibtex, ref, tags }
          bibitems.set(ref, `${bibitems.size+1}`)
        } else {
          bibitem_asts[ref].tags = array.unique([...bibitem_asts[ref].tags, ...tags])
        }
        ast.push({ type: 'text', text: text.substring(i, m.index), tags })
        ast.push({ type: 'cite', text: `[${bibitems.get(ref) as string}]`, ref, tags })
        i = m.index + m[0].length
      } else if (reftype === 'figref') {
        if (!(ref in figure_asts)) {
          figure_asts[ref] = { type: 'figure', text: '', ref, tags }
        } else {
          figure_asts[ref].tags = array.unique([...figure_asts[ref].tags, ...tags])
        }
        ast.push({ type: 'text', text: text.substring(i, m.index), tags })
        ast.push({ type: 'figref', text: `Fig.`, ref, tags })
        i = m.index + m[0].length
      }
    }
    ast.push({ type: 'text', text: `${text.substring(i)}`, tags })
    if (tags.includes('abstract')) ast[ast.length-1].text += ' '
    else if (tags.includes('introduction')) ast[ast.length-1].text += `\n\n`
    else if (tags.includes('methods')) ast[ast.length-1].text += `\n\n`
  }

  ast.push(...dict.values(figure_asts))
  ast.push(...dict.values(bibitem_asts))

  return { ast, bibitems, figures }
}
