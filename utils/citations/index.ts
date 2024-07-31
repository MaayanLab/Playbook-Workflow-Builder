import * as array from '@/utils/array'
import * as dict from '@/utils/dict'
import doiMappings from './citations'

function convertDOI(text: string) {
  const m = /^doi:(.+)$/.exec(text)
  if (m !== null && m[1] in doiMappings) {
    return `${doiMappings[m[1] as keyof typeof doiMappings]} ${text}`
  } else if (m !== null) {
    console.warn(`${m[1]} not registered, call 'npm run codegen:cite -- ${m[1]}'`)
  }
  return text
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
 */
export default function extractCitations(texts: { text?: string, tags: string[] }[]) {
  const ast: (
    { type: 'text', text: string, tags: string[] }
    | { type: 'cite', text: string, ref: string, tags: string[] }
    | { type: 'bibitem', text: string, ref: string, tags: string[] }
    | { type: 'ref', text: string, ref: string, tags: string[] }
    | { type: 'figref', text: string, ref: string, tags: string[] }
    | { type: 'figure', text: string, ref: string, tags: string[] }
  )[] = []
  const figures = new Map<string, string>()
  const figure_asts: Record<string, { type: 'figure', ref: string, text: string, tags: string[] }> = {}
  const bibitems = new Map<string, string>()
  const bibitem_asts: Record<string, { type: 'bibitem', ref: string, text: string, tags: string[] }> = {}
  for (const { text, tags } of texts) {
    if (!text) continue
    if (tags.includes('legend')) {
      const ref = tags.join('-')
      if (!figures.has(ref)) {
        figure_asts[ref] = { type: 'figure', text: `Figure ${figures.size+1}`, ref, tags }
        ast.push(figure_asts[ref])
        figures.set(ref, `${figures.size+1}`)
      } else {
        figure_asts[ref].tags = array.unique([...figure_asts[ref].tags, ...tags])
      }
    }
    const expr = /\\(figref|ref)\{(.+?)\}/g
    let i = 0
    let m
    while ((m = expr.exec(text)) !== null) {
      const reftype = m[1]
      const currentCitation = m[2]
      if (reftype === 'ref') {
        if (!bibitems.has(currentCitation)) {
          bibitem_asts[currentCitation] = { type: 'bibitem', text: `${bibitems.size+1}. ${convertDOI(currentCitation)}`, ref: currentCitation, tags }
          bibitems.set(currentCitation, `${bibitems.size+1}`)
        } else {
          bibitem_asts[currentCitation].tags = array.unique([...bibitem_asts[currentCitation].tags, ...tags])
        }
        ast.push({ type: 'text', text: text.substring(i, m.index), tags })
        ast.push({ type: 'cite', text: `[${bibitems.get(currentCitation) as string}]`, ref: currentCitation, tags })
        i = m.index + m[0].length
      } else if (reftype === 'figref') {
        if (!bibitems.has(currentCitation)) {
          figure_asts[currentCitation] = { type: 'figure', text: `${figures.size+1}. ${convertDOI(currentCitation)}`, ref: currentCitation, tags }
          figures.set(currentCitation, `${figures.size+1}`)
        } else {
          figure_asts[currentCitation].tags = array.unique([...figure_asts[currentCitation].tags, ...tags])
        }
        ast.push({ type: 'text', text: text.substring(i, m.index), tags })
        ast.push({ type: 'figref', text: `Fig. ${figures.get(currentCitation) as string}`, ref: currentCitation, tags })
        i = m.index + m[0].length
      }
    }
    ast.push({ type: 'text', text: `${text.substring(i)}`, tags })
    if (tags.includes('abstract')) ast.push({ type: 'text', text: ` `, tags })
    else if (tags.includes('introduction')) ast.push({ type: 'text', text: `\n\n`, tags })
    else if (tags.includes('methods')) ast.push({ type: 'text', text: `\n\n`, tags })
  }
  ast.push(...dict.values(bibitem_asts))
  return { ast, bibitems, figures }
}
