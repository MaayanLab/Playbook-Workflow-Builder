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
 * Using some thing [\\ref{Some Citation},\\ref{Some Other Citation}]... And also [\\ref{Some Citation}].
 * ```
 * Output:
 * ```
 * Using some thing [1,2]... And also [1].
 *
 * 1. Some Citation
 * 2. Some Other Citation
 * ```
 */
export default function extractCitations(texts: Record<string, Record<string, string | undefined>>) {
  const citations = new Map<string, string>()
  const ast: Record<string, ({ id: string, type: 'text', text: string } | { id: string, type: 'cite', ref: string } | { type: 'ref', ref: string, text: string })[]> = {}
  for (const id in texts) {
    const T = texts[id]
    for (const section in T) {
      const text = T[section as keyof typeof T]
      if (!text) continue
      if (!(section in ast)) ast[section] = []
      const expr = /\\ref\{(.+?)\}/g
      let i = 0
      let m
      while ((m = expr.exec(text)) !== null) {
        const currentCitation = m[1]
        if (!citations.has(currentCitation)) {
          citations.set(currentCitation, `${citations.size}`)
        }
        ast[section].push({ id, type: 'text', text: text.substring(i, m.index) })
        ast[section].push({ id, type: 'cite', ref: citations.get(currentCitation) as string })
        i = m.index + m[0].length
      }
      ast[section].push({ id, type: 'text', text: text.substring(i) })
    }
  }
  if (citations.size > 0) {
    ast['references'] = []
    for (const [currentCitationNum, currentCitation] of citations.entries()) {
      ast['references'].push({ type: 'ref', text: convertDOI(currentCitation), ref: currentCitationNum})
    }
  }
  return ast
}
