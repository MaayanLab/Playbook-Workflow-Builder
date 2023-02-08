import styles from './breadcrumbs.module.css'
import type { Icon as IconT } from '@/icons'
import Icon from '@/app/components/icon'
import * as dict from '@/utils/dict'
import { ensureArray } from '@/utils/array'

export default function Breadcrumbs(
  { graph, onclick: _onclick }: {
    graph: Array<{
      id: string,
      kind: 'data' | 'process',
      label: string,
      color: string,
      parents?: string[],
      icon?: IconT,
      content?: string
    }>,
    onclick?: (evt: React.MouseEvent, node: string) => void
  }
) {
  const onclick = _onclick === undefined ? () => {} : _onclick
  const g = dict.init(
    graph.map((value, index) => ({
      key: value.id,
      value: {
        ...value,
        parents: (value.parents || []),
        index,
        x: 2 * (index + 1) - 1,
      },
    }))
  )
  const e: Array<{ src: string, dst: string, dist: number, side: number }> = []
  let side = -1 // alternate bezier curve position (top/bottom)
  for (const i in g) {
    const { id, index, parents } = g[i]
    // edges
    for (const parent of parents) {
      if (!(id in g) || !(parent in g)) {
        console.error({ id: id, parent, g, graph })
        continue
      }
      // how many hops is this node away from the parent its connected to
      const dist = Math.abs(index - g[parent].index)
      // add the edge from the parent to the node with the distance and side
      e.push({
        src: parent,
        dst: id,
        dist,
        side,
      })
      // if the distance is > 1 hop (bezier instead of line), we'll alternate the side
      if (dist > 1) side *= -1
    }
  }
  
  let w = 2 * graph.length
  let h = 0.6 * Math.max(6, graph.length)
  return (
    <svg
      className="flex-grow"
      viewBox={`0 ${-h/2} ${w} ${h}`}
      preserveAspectRatio="xMinYMid meet"
    >
      {e.map((d) => {
        const { src, dst, side, dist } = d
        let path
        if (dist === 0) { // self
          // self loop
          throw new Error('NotImplemented')
        } else if (dist === 1) { // direct neighbors
          // draw line
          path = `M${g[src].x},0 L${g[dst].x},0`
        } else {
          // draw bezier curve
          path = `M${g[src].x},0 C${g[src].x},${0.4*side*dist} ${g[dst].x},${0.4*side*dist} ${g[dst].x},0`
        }
        return (
          <path
            key={`${src}__${dst}:${side}`}
            className={styles.edges}
            stroke="black"
            strokeWidth={0.05}
            fill="none"
            d={path}
          />
        )
      })}
      {dict.values(g).map((d) => {
        const { id, label, x, color, kind } = d
        const title = label || ensureArray('icon' in d ? d.icon : []).map(({ title }) => title).join(': ') || id
        return (
          <g
            key={`${id}`}
            className={`${styles.nodes} ${styles.data}`}
            transform={`translate(${x} 0)`}
            onClick={(evt) => onclick(evt, d.id)}
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
            {'icon' in d ? (
              <g transform={`scale(0.035 0.035) translate(-12 -12)`}>
                <Icon icon={d.icon} without_svg />
              </g>
            ) : (
              <text
                x={0}
                y={0}
                fontSize="0.5px"
                dominantBaseline="middle"
                textAnchor="middle"
              >{d.content}</text>
            )}
            <title>{title}</title>
          </g>
        )
      })}
    </svg>
  )
}
