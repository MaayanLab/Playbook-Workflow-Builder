import krg from '@/app/krg'
import { useRouter } from 'next/router'
import React from 'react'
import Icon from '@/app/components/icon'
import * as d3 from 'd3'
import * as array from '@/utils/array'
import { func_icon, start_icon, variable_icon } from '@/icons'

function Graph<
  N extends { id: string, x?: number, y?: number },
  L extends { id: string, source: N, target: N }
>({ origin, nodes, links, onClick }: { origin: string, nodes: N[], links: L[], onClick?: (node: N, evt: any) => void }) {
  const svgRef = React.useRef<any>(null)
  const [width, height] = [800, 600]

  React.useEffect(() => {
    if (!svgRef.current) return

    const N = [...nodes]
    const L = [...links]

    // capture origin node based on id
    let originNode = { id: origin } as N
    // recover existing node attributes after state updates
    for (const n of N) {
      if (n.id === origin) originNode = n
      const n_ref = svgRef.current.getElementById(`explore-${n.id}`)
      if (!n_ref) continue
      const transform = n_ref.getAttribute('transform')
      if (!transform) continue
      const m = /translate\((?<x>\d+(\.\d+)?),?\s+(?<y>\d+(\.\d+)?)\)/.exec(transform)
      if (!m || !m.groups) continue
      Object.assign(n, { x: Number(m.groups.x), y: Number(m.groups.y) })
    }
    // origin node initializes in the middle if not yet initialized
    if (originNode.x === undefined || originNode.y === undefined) {
      Object.assign(originNode, { x: width / 2, y: height / 2 })
    }
    // other nodes initialize near the origin (node that was clicked
    for (const n of N) {
      if (n.x === undefined || n.y === undefined) {
        Object.assign(n, {
          x: Number(originNode.x) + 50*(Math.random()-0.5),
          y: Number(originNode.y) + 50*(Math.random()-0.5),
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
  }, [svgRef.current, nodes, links, origin])

  return (
    <svg
      ref={svgRef}
      className="w-full"
      viewBox={`0 0 ${width} ${height}`}
      preserveAspectRatio="xMidYMid"
    >
      <g>
        {links.map(link => (
          <line
            key={`${link.source.id}__${link.target.id}`}
            id={`explore-${link.source.id}__${link.target.id}`}
            stroke="#999"
            strokeOpacity={0.6}
            strokeWidth={1}
          />
        ))}
        {nodes.map(node => (
          <g
            key={node.id}
            id={`explore-${node.id}`}
            onClick={evt => {
              if (onClick !== undefined) {
                onClick(node, evt)
              }
            }}
          >
            {node.id === 'Start' ? (
              <>
                <circle r="14" fill="#ddd" />
                <g transform='translate(-12, -12)'>
                  <Icon
                    icon={start_icon}
                    without_svg
                  />
                </g>
              </>
            ) : null}
            {krg.getDataNode(node.id) !== undefined ? (
              <>
                <circle r="14" fill={krg.getDataNode(node.id).meta.color || '#ddd'} />
                <g transform='translate(-12, -12)'>
                  <Icon
                    icon={krg.getDataNode(node.id).meta.icon || variable_icon}
                    title={krg.getDataNode(node.id).meta.label}
                    without_svg
                  />
                </g>
              </>
            ) : null}
            {krg.getProcessNode(node.id) !== undefined ? (
              <g transform='scale(0.5, 0.5)'>
                <rect transform='translate(-15, -15)' width="30" height="30" fill={krg.getProcessNode(node.id).meta.color || '#ddd'} />
                <g transform='translate(-12, -12)'>
                  <Icon
                    icon={krg.getProcessNode(node.id).meta.icon || func_icon}
                    title={krg.getProcessNode(node.id).meta.label}
                    without_svg
                  />
                </g>
              </g>
            ) : null}
          </g>
        ))}
      </g>
    </svg>
  )
}

export default function Explore() {
  const router = useRouter()
  const path = [...array.ensureArray(router.query.all)]
  if (path.length === 0) path.push('Start')
  const links_ = [
    ...path.slice(1).map((p, ind) => ({ source: path[ind], target: p })),
    ...path
      .flatMap(node => krg.getNextProcess(node === 'Start' ? '' : node))
      .flatMap(edge => [
        ...Object.values(edge.inputs).length === 0 ?
          [{ source: 'Start', target: edge.spec }]
        : Object.values(edge.inputs).map(input => 
            ({ source: input.spec, target: edge.spec })
          ),
        { source: edge.spec, target: edge.output.spec },
      ])
  ].reduce((agg, link) => ({ ...agg, [`${link.source}__${link.target}`]: link }), {} as Record<string, { source: string, target: string }>)
  const nodes_ = Object.values(links_).reduce(
    (agg, { source, target }) => ({
      ...agg,
      [source]: { ...(agg[source] || {}), id: source },
      [target]: { ...(agg[target] || {}), id: target },
    }),
    {} as Record<string, { id: string }>
  )
  const nodes = Object.values(nodes_)
  const links = Object.keys(links_).map((id) => ({ ...links_[id], id, source: nodes_[links_[id].source], target: nodes_[links_[id].target] }))
  return (
    <div className="col">
      <Graph
        origin={path[path.length-1]}
        nodes={nodes}
        links={links}
        onClick={(node, evt) => {
          router.push({
            pathname: '/explore/[...all]',
            query: {
              all: [...path, node.id]
            },
          }, `/explore/${[...path, node.id].join('/')}`, { shallow: true })
        }}
      />
    </div>
  )
}