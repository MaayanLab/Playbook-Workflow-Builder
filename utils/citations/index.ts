import doiMappings from './citations'
import * as dict from '@/utils/dict'

type CitableT = (
  | { type: 'section', section: string, parts: CitableT[] }
  | { type: 'text', value: string }
  | { type: 'doi', value: string }
  | { type: 'cite', value: string }
  | { type: 'tableRef', value: string }
  | { type: 'figureRef', value: string }
)

export default class Citable {
  constructor(
    public parts: CitableT[] = [],
  ) {}

  static text(strings: TemplateStringsArray, ...variables: (Citable | string)[]) {
    return new Citable(strings.flatMap((value, i) => {
      if (i > 0) {
        const variable = variables[i-1]
        if (variable instanceof Citable) {
          return [...variable.parts, { type: 'text', value }]
        } else {
          return [{ type: 'text', value: variable }, { type: 'text', value }]
        }
      } else {
        return [{ type: 'text', value }]
      }
    }))
  }
  
  static ensure(value: Citable | string): Citable {
    if (value instanceof Citable) return value
    else return new Citable([{ type: 'text', value }])
  }
  static concat(texts: (Citable | string | undefined)[]) {
    return new Citable(
      texts.flatMap((text) => text ? Citable.ensure(text).parts : [])
    )
  }
  isEmpty(): boolean {
    return this.parts.length === 0
  }
  static empty(): Citable {
    return new Citable()
  }
  static tableRef(value: string): Citable {
    return new Citable([{ type: 'tableRef', value }])
  }
  static figureRef(value: string): Citable {
    return new Citable([{ type: 'figureRef', value }])
  }
  static cite(value: string): Citable {
    return new Citable([{ type: 'cite', value }])
  }
  static doi(value: string): Citable {
    return new Citable([{ type: 'doi', value }])
  }
  section(section: string) {
    return new Citable([{ type: 'section', section, parts: this.parts }])
  }
  text(section?: string, format = 'plain' as 'plain' | 'latex'): string {
    let citations = new Map()
    let figures = new Map()
    let tables = new Map()
    const Q = [...this.parts.map(citable => ['' as string, citable] as const)].reverse()
    const S: Record<string, string> = {['undefined']: ''}
    while (Q.length > 0) {
      const [currentSection, part] = Q.pop() as [string, CitableT]
      if (part.type === 'section') {
        for (const p of part.parts.reverse())
          Q.push([part.section, p] as const)
      } else if (part.type === 'text') {
        if (section === undefined || currentSection === section) {
          if (!(currentSection in S)) S[currentSection] = ''
          S[currentSection] += part.value
        }
      } else if (part.type === 'cite' || part.type === 'doi') {
        const value = part.type === 'doi' ? `doi:${part.value}` : part.value
        if (!citations.get(value)) citations.set(value, citations.size)
        if (section === undefined || currentSection === section) {
          if (!(currentSection in S)) S[currentSection] = ''
          S[currentSection] += format === 'latex' ? `\\cite{${citations.get(value)}}` : `${citations.get(value)}`
          if (!('references' in S)) S['references'] = ''
          const inner = part.type === 'doi' ? doiMappings[part.value as keyof typeof doiMappings] : part.value
          S['references'] += format === 'latex' ? `\\bibitem{${inner}}\n` : `[${inner}]\n`
        }
      } else if (part.type === 'figureRef') {
        if (!figures.get(part.value)) figures.set(part.value, figures.size)
        if (section === undefined || currentSection === section) {
          if (!(currentSection in S)) S[currentSection] = ''
          S[currentSection] += format === 'latex' ? `\\ref{figure-${figures.get(part.value)}}` : `Fig. ${figures.get(part.value)}`
        }
      }
      else if (part.type === 'tableRef') {
        if (!tables.get(part.value)) tables.set(part.value, tables.size)
        if (section === undefined || currentSection === section) {
          if (!(currentSection in S)) S[currentSection] = ''
          S[currentSection] += format === 'latex' ? `\\ref{table-${tables.get(part.value)}}` : `Table ${tables.get(part.value)}`
        }
      }
    }
    if (!section) {
      return dict.items(S).flatMap(({ key, value }) => key ? `# ${key}\n${value}` : value).join('\n\n')
    } else {
      return S[section]
    }
  }
}
