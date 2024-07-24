import React from 'react'
import extractCitations from "@/utils/citations";
import dynamic from 'next/dynamic';
import { references_icon } from '@/icons';
import { useStory } from '../story';

const Linkify = dynamic(() => import('@/utils/linkify'))
const Icon = dynamic(() => import('@/app/components/icon'))

export function abstractText(story: ReturnType<typeof extractCitations>) {
  return story.ast.filter(part => part.tags.includes('abstract')).map(part => part.text).join('')
}

export function Abstract({ story }: { story: ReturnType<typeof extractCitations> }) { 
  const storyFiltered = React.useMemo(() => story.ast.filter((part, i) => part.tags.includes('abstract')), [story.ast])
  if (!storyFiltered.length) return null
  return (
    <p className="prose-lg text-justify mt-1">
      {storyFiltered.map((part, i) =>
        part.type === 'text' ? <Linkify key={i}>{part.text}</Linkify>
        : part.type === 'cite' ? <span key={i}> [<a key={i} href={`#${part.ref}`}>{story.bibitems.get(part.ref)}</a>]</span>
        : null
      )}
    </p>
  )
}

export function AbstractPart({ id, story }: { id: string, story: ReturnType<typeof extractCitations> }) { 
  const storyFiltered = React.useMemo(() => story.ast.filter((part, i) => part.tags.includes('abstract') && part.tags.includes(id)), [story.ast])
  if (!storyFiltered.length) return null
  return (
    <p className="prose max-w-none text-justify">
      {storyFiltered.map((part, i) =>
        part.type === 'text' ? <Linkify key={i}>{part.text}</Linkify>
        : part.type === 'cite' ? <span key={i}> [<a key={i} href={`#${part.ref}`}>{story.bibitems.get(part.ref)}</a>]</span>
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
        : part.type === 'cite' ? <span key={i}> [<a key={i} href={`#${part.ref}`}>{story.bibitems.get(part.ref)}</a>]</span>
        : null
      )}
    </p>
  )
}

export function FigureCaption({ id, story }: { id: string, story: ReturnType<typeof extractCitations> }) { 
  const storyFiltered = React.useMemo(() => story.ast.filter((part, i) => part.tags.includes('legend') && part.tags.includes(id)), [id, story.ast])
  if (!storyFiltered.length) return null
  return (
    <div className="prose max-w-none">
      {storyFiltered.map((part, i) =>
        part.type === 'figure' ? <strong key={i}>{part.text}.&nbsp;</strong>
        : part.type === 'text' ? <Linkify key={i}>{part.text}</Linkify>
        : part.type === 'cite' ? <span key={i}> [<a key={i} href={`#${part.ref}`}>{story.bibitems.get(part.ref)}</a>]</span>
        : null
      )}
    </div>
  )
}

export function References() {
  const story = useStory()
  return  (
    <div className="flex-grow flex-shrink bp5-card p-0">
      <div className="p-3">
        <div className="flex flex-row gap-2">
          <Icon icon={references_icon} className="fill-black dark:fill-white" />
          <h2 className="bp5-heading w-full">
            References
          </h2>
        </div>
        <div className="prose max-w-none">
          <p className="text-sm text-justify whitespace-pre-line my-2 flex flex-col">{story.ast.map((part, i) =>
            part.type === 'bibitem' ? <a key={i} id={part.ref} className="no-underline text-neutral-focus hover:text-neutral hover:no-underline hover:cursor-text"><Linkify>{part.text}</Linkify></a>
              : null
          )}</p>
        </div>
      </div>
    </div>
  )
}
