import React from 'react'
import dynamic from 'next/dynamic'
import extractCitations from '@/utils/citations'
import type { DataMetaNode, ProcessMetaNode, Story } from '@/spec/metanode'
import * as dict from '@/utils/dict'

const Linkify = dynamic(() => import('@/utils/linkify'))
const Bibliography = dynamic(() => import('@/app/fragments/bibliography'))

const sections: { key: keyof Story, label: string }[] = [
  { key: 'abstract', label: 'Abstract' },
  { key: 'introduction', label: 'Introduction' },
  { key: 'methods', label: 'Methods' },
  { key: 'results', label: 'Results' },
  { key: 'figureLegend', label: 'Figure legend' },
  { key: 'legend', label: 'Figure legend' },
  { key: 'tableLegend', label: 'Table legend' },
]

/**
 * The story of a component as it'd be written in a report, "baked" without any data (like the
 *  table export, see cli/build-table): what depends on the data reads as the story's defaults.
 *  The references to figures (of the inputs & of this step's output) are placeholders.
 */
function bakeStory(metanode: ProcessMetaNode) {
  if (!metanode.story) return undefined
  const input_refs = dict.init(dict.items(metanode.inputs ?? {}).map(({ key, value }) => {
    const ref = `\\figref{input:${String(key)}}`
    return { key, value: Array.isArray(value) ? [ref] : ref }
  }))
  try {
    return metanode.story({ input_refs, output_ref: '\\figref{output}' } as any) as Story
  } catch (e) {
    console.warn(`${metanode.spec}: its story failed without data`, e)
    return undefined
  }
}

function figureLabel(ref: string) {
  if (ref === 'output') return `this step's figure`
  const m = /^input:(.+)$/.exec(ref)
  return m ? <>the <em>{m[1]}</em> input's figure</> : ref
}

export default function ComponentStory({ metanode }: { metanode: DataMetaNode | ProcessMetaNode }) {
  const story = React.useMemo(() => {
    if (metanode.kind !== 'process') return undefined
    const baked = bakeStory(metanode)
    if (!baked) return undefined
    const present = sections.filter(({ key }) => typeof baked[key] === 'string' && (baked[key] as string).trim() !== '')
    if (present.length === 0) return undefined
    const citations = extractCitations(present.map(({ key }) => ({ text: baked[key], tags: [key] })))
    // the index of each section's last text
    const lastText: Partial<Record<keyof Story, number>> = {}
    citations.ast.forEach((part, i) => { if (part.type === 'text') for (const tag of part.tags) lastText[tag as keyof Story] = i })
    return { present, lastText, ...citations }
  }, [metanode])
  if (!story) return null
  const bibitems = story.ast.flatMap(part => part.type === 'bibitem' ? [part] : [])
  return (
    <div className="flex flex-col gap-3 my-2">
      <h2 className="prose max-w-none text-xl font-bold">Story</h2>
      <p className="prose max-w-none text-sm opacity-70 -mt-2">How this step is described in a report (before it has any data: what's specific to the data reads as the defaults, and figures are placeholders).</p>
      {story.present.map(({ key, label }) => (
        <div key={key} className="prose max-w-none">
          <strong>{label}</strong>:{' '}
          <span className="whitespace-pre-line">
            {story.ast.map((part, i) =>
              !part.tags.includes(key) ? null
              // (extractCitations separates a report's paragraphs, each section stands alone here)
              : part.type === 'text' ? <Linkify key={i}>{i === story.lastText[key] ? part.text.trimEnd() : part.text}</Linkify>
              : part.type === 'cite' ? <span key={i}>[<a href={`#${part.ref}`}>{story.bibitems.get(part.ref)}</a>]</span>
              : part.type === 'figref' ? <span key={i} className="badge badge-ghost align-baseline">{figureLabel(part.ref)}</span>
              : null
            )}
          </span>
        </div>
      ))}
      {bibitems.length ? (
        <div className="prose max-w-none">
          <strong>References</strong>
          <Bibliography className="text-sm my-2" items={bibitems} />
        </div>
      ) : null}
    </div>
  )
}
