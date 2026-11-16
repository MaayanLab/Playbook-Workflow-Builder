import React from 'react'
import { ExLink } from '@/app/fragments/ex-router'
import dynamic from 'next/dynamic'
import { useVirtualizer } from '@tanstack/react-virtual'
import krg from "@/app/krg"
import { func_icon, variable_icon } from '@/icons'
import { MetaNode, DataMetaNode, ProcessMetaNode } from '@/spec/metanode'
import { metanodeType as getMetaNodeType } from '@/app/metanode-type'

const Markdown = dynamic(() => import('@/app/components/Markdown'))
const Icon = dynamic(() => import('@/app/components/icon'))

export function metanodeType(metanode: MetaNode) {
  return getMetaNodeType(krg, metanode)
}

type SortKey = 'type' | 'label'

/** The element which scrolls `el` (the shell's page, see @/app/fragments/playbook/layout) */
function scrollParent(el: HTMLElement | null) {
  for (let parent = el?.parentElement; parent; parent = parent.parentElement) {
    const { overflowY } = getComputedStyle(parent)
    if (overflowY === 'auto' || overflowY === 'scroll') return parent
  }
  return document.scrollingElement as HTMLElement | null
}

const compare = (a: string, b: string) => a.localeCompare(b, undefined, { sensitivity: 'base', numeric: true })

// icon | type | label (| description, from md up: below the type & label on small screens)
const columns = 'grid grid-cols-[2.5rem_7rem_minmax(0,1fr)] md:grid-cols-[3rem_11rem_minmax(0,1fr)_minmax(0,2fr)] gap-x-3'

/**
 * A searchable, sortable listing of components: its toolbar (search & column headers) stays
 *  visible while scrolling, and only the rows near the visible part of the page are rendered.
 */
export default function MetaNodeListing({ metanodes }: { metanodes: { metanode: (DataMetaNode | ProcessMetaNode), multi?: boolean, type?: string, }[] }) {
  const [search, setSearch] = React.useState('')
  const [sort, setSort] = React.useState<{ key: SortKey, desc: boolean }>({ key: 'type', desc: false })
  const rows = React.useMemo(() => metanodes.map(({ metanode, type }) => {
    const rowType = type ?? metanodeType(metanode)
    return {
      key: `${rowType}-${metanode.spec}`,
      metanode,
      type: rowType,
      searchable: [rowType, metanode.meta.label, metanode.spec, metanode.meta.description ?? ''].join(' ').toLowerCase(),
    }
  }), [metanodes])
  const shown = React.useMemo(() => {
    const terms = search.toLowerCase().split(/\s+/).filter(term => term !== '')
    const value = (row: typeof rows[number], key: SortKey) => key === 'type' ? row.type : row.metanode.meta.label
    const other: SortKey = sort.key === 'type' ? 'label' : 'type'
    return rows
      .filter(row => terms.every(term => row.searchable.includes(term)))
      .sort((a, b) => (sort.desc ? -1 : 1) * (compare(value(a, sort.key), value(b, sort.key)) || compare(value(a, other), value(b, other))))
  }, [rows, search, sort])

  // the rows are positioned within the page's scroll container
  const listRef = React.useRef<HTMLDivElement>(null)
  const [scroller, setScroller] = React.useState<HTMLElement | null>(null)
  const [scrollMargin, setScrollMargin] = React.useState(0)
  React.useLayoutEffect(() => { setScroller(scrollParent(listRef.current)) }, [])
  React.useLayoutEffect(() => {
    if (!scroller) return
    // where the list starts in the scroll container (what's above it can change size, e.g. alerts)
    const update = () => {
      if (!listRef.current) return
      setScrollMargin(listRef.current.getBoundingClientRect().top - scroller.getBoundingClientRect().top + scroller.scrollTop)
    }
    update()
    const observer = new ResizeObserver(update)
    for (const child of scroller.children) observer.observe(child)
    return () => observer.disconnect()
  }, [scroller])
  const virtualizer = useVirtualizer({
    count: shown.length,
    getScrollElement: () => scroller,
    estimateSize: () => 64,
    overscan: 8,
    scrollMargin,
    getItemKey: index => shown[index].key,
  })

  const SortHeader = ({ sortKey, children }: React.PropsWithChildren<{ sortKey: SortKey }>) => (
    <button
      className="flex flex-row items-center gap-1 font-bold text-left hover:underline"
      aria-sort={sort.key === sortKey ? (sort.desc ? 'descending' : 'ascending') : 'none'}
      onClick={() => setSort(current => current.key === sortKey ? { key: sortKey, desc: !current.desc } : { key: sortKey, desc: false })}
    >
      {children}
      <span className={sort.key === sortKey ? '' : 'opacity-30'}>{sort.key === sortKey && sort.desc ? '▼' : '▲'}</span>
    </button>
  )

  return (
    <div className="w-full flex flex-col">
      <div className="sticky top-0 z-10 bg-base-100 border-b border-current/20 flex flex-col gap-2 pt-2 pb-1">
        <div className="flex flex-row items-center gap-3">
          <input
            type="search"
            className="input input-bordered input-sm grow"
            placeholder="Search components..."
            value={search}
            onChange={evt => setSearch(evt.currentTarget.value)}
          />
          <span className="whitespace-nowrap text-sm opacity-70">{shown.length === rows.length ? rows.length : `${shown.length} of ${rows.length}`}</span>
        </div>
        <div className={`${columns} px-2 text-sm`} role="row">
          <span>Icon</span>
          <SortHeader sortKey="type">Type</SortHeader>
          <SortHeader sortKey="label">Label</SortHeader>
          <span className="hidden md:block font-bold">Description</span>
        </div>
      </div>
      {shown.length === 0 ? <div className="prose max-w-none p-4">No components match your search.</div> : null}
      <div ref={listRef} className="relative w-full" style={{ height: virtualizer.getTotalSize() }}>
        {virtualizer.getVirtualItems().map(item => {
          const { metanode, type } = shown[item.index]
          return (
            <div
              key={item.key}
              data-index={item.index}
              ref={virtualizer.measureElement}
              className={`${columns} absolute left-0 top-0 w-full items-start px-2 py-3 border-b border-current/10 hover:bg-current/5`}
              style={{ transform: `translateY(${item.start - virtualizer.options.scrollMargin}px)` }}
            >
              <Icon icon={metanode.meta.icon ? metanode.meta.icon : metanode.kind === 'data' ? variable_icon : func_icon} />
              <span className="text-sm">{type}</span>
              <ExLink
                shallow
                href={`/components/${encodeURIComponent(metanode.spec)}`}
                className="link prose max-w-none font-bold whitespace-pre-wrap [overflow-wrap:anywhere]"
              >{metanode.meta.label}</ExLink>
              <div className="col-start-2 col-span-2 md:col-start-auto md:col-span-1 prose prose-sm max-w-none whitespace-pre-wrap [overflow-wrap:anywhere]">
                <Markdown components={{
                  p(props: React.PropsWithChildren<{}>) {
                    return <span {...props}>{props.children}</span>
                  }
                }}>{metanode.meta.description}</Markdown>
              </div>
            </div>
          )
        })}
      </div>
    </div>
  )
}

export function MetaNodeCard({ metanode }: { metanode: DataMetaNode | ProcessMetaNode }) {
  return (
    <div className="flex flex-row items-center gap-2 m-2">
      <Icon icon={metanode.meta.icon ? metanode.meta.icon : metanode.kind === 'data' ? variable_icon : func_icon} size={3} title={null} />
      <div className="flex flex-col prose max-w-none">
        <h1 className="m-0 whitespace-pre-wrap">{metanode.meta.label}</h1>
        <h2 className="m-0 text-gray-500">{metanode.kind === 'data' ? 'Data' : 'prompt' in metanode ? 'Prompt' : 'Resolver'} Node</h2>
      </div>
    </div>
  )
}
