import React from 'react'
import { Map } from 'immutable'
import * as dict from '@/utils/dict'

export type BreadcrumbNode = {
  id: string,
  index: number,
  parents?: string[],
  onClick?: React.MouseEventHandler<SVGGElement>,
}

/**
 * Breadcrumb layout algorithm -- at a high level, it's a depth first traversal
 *  through node children where steps down the tree expand in the x direction
 *  and steps over in the tree expand in the y direction.
 */
function layout(graph: Map<string, BreadcrumbNode>) {
  // Construct C -- a mapping from id to children
  //         & G -- a mapping from id to node
  const C: Record<string, string[]> = {}
  if (graph.size === 0) return { P: {}, E: [], W: 0, H: 0 }
  const G = dict.init(
    graph.valueSeq().sort((a, b) => a.index - b.index).toArray()
      .map((node) => {
        (node.parents??[]).forEach(pid => {
          if (!(pid in C)) C[pid] = []
          C[pid].push(node.id)
        })
        return { key: node.id, value: node }
      })
  )
  // DFS Through C to create B -- the ordered breadcrumbs with x/y postiions
  //                       & E -- the edges between breadcrumbs
  // W & H keep track of the width/height of the layout
  let W = 1, H = 0
  const Q = [{id: 'start', parent: undefined as string | undefined, x: W-1, y: H++}]
  const P: Record<string, {x: number, y: number}> = {}
  const E: {id: string, src: string, dst: string}[] = []
  while (Q.length > 0) {
    const el = Q.pop()
    if (!el) break
    if (el.y === -1) {
      if (el.id in P) el.y = P[el.id].y
      else if ((G[el.id].parents||[]).some(pid => !(pid in P))) {
        // defer until all parents have been shown
        Q.splice(Q.length-2, 0, el)
        continue
      } else el.y = H++
    }
    if (el.parent) {
      E.push({
        id: `${el.parent}__${el.id}__${el.x},${el.y}`,
        src: el.parent,
        dst: el.id,
      })
    }
    W = Math.max(el.x+1, W)
    if (el.id in P) continue
    P[el.id] = { x: el.x, y: el.y }
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
  return { P, E, W, H }
}

const BreadcrumbContext = React.createContext(undefined as {
  add: (node: BreadcrumbNode) => void,
  del: (node: BreadcrumbNode) => void,
  pos: (node: BreadcrumbNode) => ({ x: number, y: number } | undefined),
} | undefined)

export function Breadcrumbs({ children }: React.PropsWithChildren<{}>) {
  const [graph, setGraph] = React.useState(Map<string, BreadcrumbNode>())
  const { P, E, top, left, w, h } = React.useMemo(() => {
    const { P, E, W, H } = layout(graph)
    let top = 1
    let left = 1
    let w = 1.5 * W - 0.5
    if (w < 0) {
      left -= w
      w = Math.abs(w)
    }
    let h = 1.5 * H + (H === 1 ? 0.5 : 0.25)
    if (h < 0) {
      top -= h
      h = Math.abs(h)
    }
    return { P, E, top, left, w, h }
  }, [graph])
  const add = React.useCallback((node: BreadcrumbNode) => {setGraph(graph => graph.set(node.id, node))}, [setGraph])
  const del = React.useCallback((node: BreadcrumbNode) => {setGraph(graph => graph.delete(node.id))}, [setGraph])
  const pos = React.useCallback((node: BreadcrumbNode) => node.id in P ? { x: P[node.id].x, y: P[node.id].y } : undefined, [P])
  return (
    <svg
      className="flex-grow"
      style={{ height: 40*h }}
      viewBox={`0 0 ${w} ${h}`}
      preserveAspectRatio="xMinYMid meet"
    >
      <g transform={`translate(${left}, ${top})`}>
        {E.map(({ id, src, dst }) => {
          return (
            <line
              key={id}
              className="stroke-black dark:stroke-white"
              strokeWidth={0.05}
              x1={P[src].x*1.25}
              y1={P[src].y*1.25}
              x2={P[dst].x*1.25}
              y2={P[dst].y*1.25}
            />
          )
        })}
        <BreadcrumbContext.Provider value={{ add, del, pos }}>
          {children}
        </BreadcrumbContext.Provider>
      </g>
    </svg>
  )
}

export function useBreadcrumbPosition(props: BreadcrumbNode) {
  const parent = React.useContext(BreadcrumbContext)
  const position = React.useMemo(() => parent?.pos(props), [parent?.pos, props])
  React.useEffect(() => {
    if (!parent) return
    parent.add(props)
    return () => {parent.del(props)}
  }, [parent?.add, parent?.del, props])
  return position
}
