import React from 'react'
import krg from '@/app/krg'
import dynamic from 'next/dynamic'
import { variable_icon } from '@/icons'
import Link from 'next/link'
import { useRouter } from 'next/router'
import { DataMetaNode } from '@/spec/metanode'

const Icon = dynamic(() => import('@/app/components/icon'))
const Card = dynamic(() => import('@blueprintjs/core').then(({ Card }) => Card))

export default function ChooseEnd() {
  const router = useRouter()
  const startNode = React.useMemo(() => {
    if (!router.isReady || !router.query || !router.query.start) return
    return krg.getDataNode(router.query.start as string)
  }, [router.query.start])
  const endNodes = React.useMemo<DataMetaNode[] | undefined>(() => {
    if (!startNode) return undefined
    const Q = krg.getNextProcess(startNode.spec)
    const D = new Set<string>()
    const C: DataMetaNode[] = []
    while (Q.length > 0) {
      const proc = Q.pop()
      if (!proc) continue
      if (D.has(proc.output.spec)) continue
      C.push(proc.output)
      D.add(proc.output.spec)
      krg.getNextProcess(proc.output.spec)
        .forEach(proc => {
          if (!D.has(proc.spec)) {
            Q.push(proc)
            D.add(proc.spec)
          }
        })
    }
    return C
  }, [startNode])
  return (
    <div className="container mx-auto flex flex-col items-center py-4">
      <h1 className="bp4-heading"><Icon icon={variable_icon} />&nbsp;End</h1>
      <div className="flex flex-row flex-wrap gap-2 justify-center">
        {startNode ? (endNodes || []).map(item => {
          return (
            <Link href={`/explore/dfs/${startNode.spec}/${item.spec}`}>
              <Card
                key={item.spec}
                className="col w-80"
                interactive
                style={{
                  backgroundColor: item.meta.color || 'lightgrey',
                }}
              >
                <h5 className="bp4-heading"><Icon title={item.meta.label} icon={item.meta.icon || variable_icon} /> {item.meta.label || ''}</h5>
                <p className="bp4-text-small">{item.meta.description || ''}</p>
              </Card>
            </Link>
          )
        }) : null}
      </div>
    </div>
  )
}