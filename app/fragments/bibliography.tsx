import React from 'react'
import dynamic from 'next/dynamic'

const Linkify = dynamic(() => import('@/utils/linkify'))

/**
 * A numbered list of references (bibitems' texts start with their number, e.g. `1. Author...`):
 *  the numbers form a column, so that each reference's lines (when it wraps) line up with its
 *  start. Left aligned, and long unbroken text (links) wraps rather than widening the page.
 */
export default function Bibliography({ items, className }: { items: { ref: string, text: string }[], className?: string }) {
  return (
    <div className={`grid grid-cols-[auto_minmax(0,1fr)] gap-x-2 gap-y-1 text-left min-w-0 ${className ?? ''}`}>
      {items.map(({ ref, text }) => {
        const m = /^(\d+\.)\s+([\s\S]*)$/.exec(text)
        return (
          <React.Fragment key={ref}>
            <span id={ref} className="text-right">{m ? m[1] : ''}</span>
            <span className="min-w-0 whitespace-pre-line [overflow-wrap:anywhere]"><Linkify>{m ? m[2] : text}</Linkify></span>
          </React.Fragment>
        )
      })}
    </div>
  )
}
