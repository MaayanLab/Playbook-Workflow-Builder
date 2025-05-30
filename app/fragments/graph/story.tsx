import React from 'react'
import extractCitations from "@/utils/citations";
import dynamic from 'next/dynamic';

const Linkify = dynamic(() => import('@/utils/linkify'))

export function Abstract({ story }: { story: ReturnType<typeof extractCitations> }) { 
  const storyFiltered = React.useMemo(() => story.ast.filter((part, i) => part.tags.includes('abstract')), [story.ast])
  if (!storyFiltered.length) return null
  return (
    <p className="prose max-w-none text-justify">
      {storyFiltered.map((part, i) =>
        part.type === 'text' ? <Linkify key={i}>{part.text}</Linkify>
        : part.type === 'cite' ? <span key={i}> [<a key={i} href={`#${part.ref}`}>{story.bibitems.get(part.ref)}</a>]</span>
        : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'figure' ? <span key={i}> <a key={i} href={`#${part.ref}`}>Fig. {story.figures.get(part.ref)?.ref}</a></span>
        : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'table' ? <span key={i}> <a key={i} href={`#${part.ref}`}>Table. {story.figures.get(part.ref)?.ref}</a></span>
        : null
      )}
    </p>
  )
}

export function Methods({ id, story }: { id: string, story: ReturnType<typeof extractCitations> }) { 
  const storyFiltered = React.useMemo(() => story.ast.filter((part, i) => part.tags.includes('methods') && part.tags.includes(id)), [id, story.ast])
  if (!storyFiltered.length) return null
  return (
    <p className="prose max-w-none text-justify">
      <strong>Method </strong>
      {storyFiltered.map((part, i) =>
        part.type === 'text' ? <Linkify key={i}>{part.text}</Linkify>
        : part.type === 'cite' ? <span key={i}> [<a href={`#${part.ref}`}>{story.bibitems.get(part.ref)}</a>]</span>
        : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'figure' ? <span key={i}> <a key={i} href={`#${part.ref}`}>Fig. {story.figures.get(part.ref)?.ref}</a></span>
        : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'table' ? <span key={i}> <a key={i} href={`#${part.ref}`}>Table. {story.figures.get(part.ref)?.ref}</a></span>
        : null
      )}
    </p>
  )
}

export function FigureCaption({ id, story }: { id: string, story: ReturnType<typeof extractCitations> }) { 
  const storyFiltered = React.useMemo(() => story.ast.filter((part, i) => (part.tags.includes('legend') || part.tags.includes('figureLegend') || part.tags.includes('tableLegend')) && part.tags.includes(id)), [id, story.ast])
  if (!storyFiltered.length) return null
  return (
    <div className="prose max-w-none">
      {story.figures.get(id)?.kind === 'figure' && <strong>Figure {story.figures.get(id)?.ref}.</strong>}
      {story.figures.get(id)?.kind === 'table' && <strong>Table {story.figures.get(id)?.ref}.</strong>}
      &nbsp;
      {storyFiltered.map((part, i) =>
        part.type === 'text' ? <Linkify key={i}>{part.text}</Linkify>
        : part.type === 'cite' ? <span key={i}> [<a href={`#${part.ref}`}>{story.bibitems.get(part.ref)}</a>]</span>
        : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'figure' ? <a key={i} href={`#${part.ref}`}>Fig. {story.figures.get(part.ref)?.ref}</a>
        : part.type === 'figref' && story.figures.get(part.ref)?.kind === 'table' ? <a key={i} href={`#${part.ref}`}>Table. {story.figures.get(part.ref)?.ref}</a>
        : null
      )}
    </div>
  )
}

export function References({ story }: { story: ReturnType<typeof extractCitations> }) {
  const storyFiltered = React.useMemo(() => story.ast.filter((part, i) => part.type === 'bibitem'), [story.ast])
  if (!storyFiltered.length) return null
  return (
    <div className="prose max-w-none">
      <strong>References</strong>
      <p className="text-sm text-justify whitespace-pre-line my-2 flex flex-col">{storyFiltered.map((part, i) =>
        part.type === 'bibitem' ? <a key={i} id={part.ref} className="no-underline text-neutral-focus hover:text-neutral hover:no-underline hover:cursor-text"><Linkify>{part.text}</Linkify></a>
        : null
      )}</p>
    </div>
  )
}
