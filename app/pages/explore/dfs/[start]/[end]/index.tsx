import React from 'react'
import krg from '@/app/krg'
import * as d3 from 'd3'
import * as dict from '@/utils/dict'
import { useRouter } from 'next/router'
import { func_icon, variable_icon } from '@/icons'
import Icon from '@/app/components/icon'

function Graph(props: {
  nodes: { id: string, x?: number, y?: number }[],
  processes: { id: string, x?: number, y?: number }[],
  links: {
    id: string,
    source: { id: string, x?: number, y?: number },
    target: { id: string, x?: number, y?: number },
  }[],
  onClick?: (node: string, evt: any) => void
}) {
  const svgRef = React.useRef<any>(null)
  const [width, height] = [800, 600]

  React.useEffect(() => {
    if (!svgRef.current) return

    const N = [...props.nodes, ...props.processes]
    const L = [...props.links]
    // recover existing node attributes after state updates
    for (const n of N) {
      const n_ref = svgRef.current.getElementById(`explore-${n.id}`)
      if (!n_ref) continue
      const transform = n_ref.getAttribute('transform')
      if (!transform) continue
      const m = /translate\((?<x>\d+(\.\d+)?),?\s+(?<y>\d+(\.\d+)?)\)/.exec(transform)
      if (!m || !m.groups) continue
      Object.assign(n, { x: Number(m.groups.x), y: Number(m.groups.y) })
    }
    // other nodes initialize near the origin (node that was clicked
    for (const n of N) {
      if (n.x === undefined || n.y === undefined) {
        Object.assign(n, {
          x: width/2 + 50*(Math.random()-0.5),
          y: height/2 + 50*(Math.random()-0.5),
        })
      }
    }

    const simulation = d3.forceSimulation(N)
      .force('charge', d3.forceManyBody().strength(-150))
      .force("link", d3.forceLink(L).distance(15))
      .force('center', d3.forceCenter(width / 2, height / 2))
      .force("x", d3.forceX())
      .force("y", d3.forceY())
      .alpha(0.1)

    simulation.on('tick', () => {
      for (const n of N) {
        const n_ref = svgRef.current.getElementById(`explore-${n.id}`)
        if (!n_ref) continue
        n_ref.setAttribute('transform', `translate(${n.x}, ${n.y})`)
      }
      for (const l of L) {
        const l_ref = svgRef.current.getElementById(`explore-${l.id}`)
        if (!l_ref) continue
        l_ref.setAttribute('x1', l.source.x)
        l_ref.setAttribute('y1', l.source.y)
        l_ref.setAttribute('x2', l.target.x)
        l_ref.setAttribute('y2', l.target.y)
      }
    })

    simulation.restart()
    // d3-controlled svg zoom
    const d3SvgRef = d3.select(svgRef.current)
    d3SvgRef
      .call(d3.zoom().on('zoom', (e) => {
        d3SvgRef
          .selectChild('g')
          .attr('transform', e.transform)
      }))

    return () => {
      simulation.stop()
    }
  }, [props.nodes, props.processes, props.links])
  return (
    <svg
      ref={svgRef}
      className="w-screen h-screen"
      viewBox={`0 0 ${width} ${height}`}
      preserveAspectRatio="xMidYMid"
    >
      <g>
        {props.links.map(link => (
          <line
            key={link.id}
            id={`explore-${link.id}`}
            stroke="#999"
            strokeOpacity={0.6}
            strokeWidth={1}
          />
        ))}
        {props.nodes.map(({ id }) => krg.getDataNode(id)).filter(node => node !== undefined).map(node => (
          <g
            key={node.spec}
            id={`explore-${node.spec}`}
            onClick={evt => {
              if (props.onClick !== undefined) {
                props.onClick(node.spec, evt)
              }
            }}
          >
            <circle r="16" fill={node.meta.color || '#ddd'} />
            <g transform='translate(-12, -12)'>
              <Icon
                icon={node.meta.icon || variable_icon}
                title={node.meta.label}
                without_svg
              />
            </g>
          </g>
        ))}
        {props.processes.map(({ id }) => krg.getProcessNode(id)).filter(node => node !== undefined).map(node => (
          <g
            key={node.spec}
            id={`explore-${node.spec}`}
            onClick={evt => {
              if (props.onClick !== undefined) {
                props.onClick(node.spec, evt)
              }
            }}
          >
            <g transform='scale(0.5, 0.5)'>
              <rect transform='translate(-16, -16)' width="32" height="32" fill={node.meta.color || '#ddd'} />
              <g transform='translate(-12, -12)'>
                <Icon
                  icon={node.meta.icon || func_icon}
                  title={node.meta.label}
                  without_svg
                />
              </g>
            </g>
          </g>
        ))}
      </g>
    </svg>
  )
}

export default function ExploreDFS() {
  const router = useRouter()
  const startNode = React.useMemo(() => {
    if (!router.isReady || !router.query || !router.query.start) return
    return krg.getDataNode(router.query.start as string)
  }, [router.query.start])
  const endNode = React.useMemo(() => {
    if (!router.isReady || !router.query || !router.query.end) return
    return krg.getDataNode(router.query.end as string)
  }, [router.query.end])
  const graph = React.useMemo(() => {
    if (!startNode || !endNode) return
    const Q = krg.getNextProcess(startNode.spec).map(proc => ({ parent: startNode, proc }))
    while (Q.length > 0) {
      const el = Q.pop()
      if (!el) continue
      const { parent, proc } = el
      if (proc.spec in B) continue
      B[proc.spec] = parent.spec
      B[proc.output.spec] = proc.spec
      krg.getNextProcess(proc.output.spec)
        .forEach(proc => Q.push({ parent: proc.output, proc }))
    }
    const nodes: Record<string, { id: string }> = {}
    const processes: Record<string, { id: string }> = {}
    const links: { id: string, source: { id: string }, target: { id: string } }[] = []
    nodes[startNode.spec] = { id: startNode.spec }
    nodes[endNode.spec] = { id: endNode.spec }
    
    let first = true
    let head = endNode.spec
    while (first || head !== startNode.spec) {
      first = false
      const proc: string = B[head]
      if (!(proc in processes)) {
        processes[proc] = { id: proc }
      } else {
        break
      }
      const out: string = B[proc]
      if (!(out in nodes)) {
        nodes[out] = { id: out }
      } else {
        break
      }
      links.push({ id: `${proc}__${out}`, source: processes[proc], target: nodes[out] })
      links.push({ id: `${head}__${proc}`, source: nodes[head], target: processes[proc] })
      head = out
    }
    return {
      nodes: dict.values(nodes),
      processes: dict.values(processes),
      links,
    }
  }, [startNode, endNode])
  
  if (!graph) return null
  return <Graph {...graph} />
}