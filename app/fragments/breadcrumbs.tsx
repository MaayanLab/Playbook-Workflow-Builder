import React from 'react'
import type { Icon as IconT } from '@/icons'
import Icon from '@/app/components/icon'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import * as d3 from 'd3'

type BreadcrumbNode = {
  id: string,
  kind: 'data' | 'process',
  label: string,
  color: string,
  parents?: string[],
  icon: IconT,
}

/**
 * Breadcrumb layout algorithm -- at a high level, it's a depth first traversal
 *  through node children where steps down the tree expand in the x direction
 *  and steps over in the tree expand in the y direction.
 */
function layout(graph: BreadcrumbNode[]) {
  // Construct C -- a mapping from id to children
  //         & G -- a mapping from id to node
  const C: Record<string, string[]> = {}
  const G = dict.init(graph.map(node => {
    (node.parents??[]).forEach(pid => {
      if (!(pid in C)) C[pid] = []
      C[pid].push(node.id)
    })
    return { key: node.id, value: node }
  }))
  // DFS Through C to create B -- the ordered breadcrumbs with x/y postiions
  //                       & E -- the edges between breadcrumbs
  // W & H keep track of the width/height of the layout
  let W = 1, H = 1
  const Q = [{id: 'start', parent: undefined as string | undefined, x: W-1, y: H-1}]
  const B: Record<string, {node: BreadcrumbNode, x: number, y: number}> = {}
  const E: {id: string, src: string, dst: string}[] = []
  while (Q.length > 0) {
    const el = Q.pop()
    if (!el) break
    if (el.y === -1) el.y = H++
    W = Math.max(el.x+1, W)
    B[el.id] = { node: G[el.id], x: el.x, y: el.y }
    if (el.parent) {
      E.push({
        id: `${el.id}__${el.id}__${el.x},${el.y}`,
        src: el.parent,
        dst: el.id,
      })
    }
    if (el.id in C) {
      // we reverse children here since the last item inserted will be processed next
      const children = [...C[el.id]]
      children.reverse()
      children.forEach((cid, rind) => {
        // maintain original children index (even though list is reversed)
        const ind = children.length - rind - 1
        Q.push({
          parent: el.id,
          id: cid,
          x: el.x+1,
          y: ind === 0 ? el.y : -1, // -1 is used to defer new branch creation
        })
      })
    }
  }
  return { B, E, W, H }
}

export default function Breadcrumbs(
  { graph, onclick: _onclick }: {
    graph: Array<BreadcrumbNode>,
    onclick?: (evt: React.MouseEvent, node: string) => void
  }
) {
  const onclick = _onclick === undefined ? () => {} : _onclick
  const { B, E, W, H } = React.useMemo(() => layout(graph), [graph])
  let w = 1.5 * Math.max(2, W)
  let h = 1.5 * Math.max(2, H)
  return (
    <svg
      className="flex-grow"
      viewBox={`-1 -1 ${w} ${h}`}
      preserveAspectRatio="xMinYMid meet"
    >
      <g>
        {E.map(({ id, src, dst }) => {
          return (
            <line
              key={id}
              stroke="black"
              strokeWidth={0.05}
              fill="none"
              x1={B[src].x*1.25}
              y1={B[src].y*1.25}
              x2={B[dst].x*1.25}
              y2={B[dst].y*1.25}
            />
          )
        })}
        {dict.values(B).map((d) => {
          const { node: { id, label, color, kind, icon }, x, y } = d
          const title = label || array.ensureArray(icon).map(({ title }) => title).join(': ') || id
          return (
            <g
              key={`${id}`}
              className="cursor-pointer"
              transform={`translate(${x*1.25} ${y*1.25})`}
              onClick={(evt) => onclick(evt, id)}
            >
              {kind === 'data' ? (
                <circle
                  fill={color}
                  stroke="black"
                  strokeWidth={0.001}
                  r={0.5}
                  cx={0}
                  cy={0}
                />
              ) : (
                <rect
                  fill={color}
                  x={-0.5}
                  y={-0.5}
                  width={1}
                  height={1}
                />
              )}
              <g transform={`scale(0.035 0.035) translate(-12 -12)`}>
                <Icon icon={icon} title={title} without_svg />
              </g>
            </g>
          )
        })}
      </g>
    </svg>
  )
}
