import useSWRImmutable from 'swr/immutable'
import { start_icon, func_icon, variable_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { useRouter } from 'next/router'
import type { Metapath } from '@/app/fragments/graph/types'
import krg from '@/app/krg'

const Breadcrumbs = dynamic(() => import('@/app/fragments/breadcrumbs'))
const Home = dynamic(() => import('@/app/fragments/playbook/home'))
const Extend = dynamic(() => import('@/app/fragments/graph/extend'))
const Suggest = dynamic(() => import('@/app/fragments/graph/suggest'))
const Cell = dynamic(() => import('@/app/fragments/graph/cell'))

export default function Graph({ graph_id, node_id, extend, suggest }: { graph_id: string, node_id: string, extend: boolean, suggest: boolean }) {
  const router = useRouter()
  const { data: metapath_, error } = useSWRImmutable<Array<Metapath>>(() => graph_id !== 'start' ? `/api/db/fpl/${graph_id}` : undefined)
  const metapath = metapath_ || []
  const head = metapath.filter(({ id }) => id === node_id)[0]
  return (
    <>
      <div className="flex w-auto h-40">
        <Breadcrumbs
          graph={[
            {
              id: 'start',
              kind: 'data',
              label: 'Start',
              color: node_id === 'start' ? '#B3CFFF' : 'lightgrey',
              icon: [start_icon],
              parents: [],
            },
            ...metapath.flatMap((head, i) => {
              const process = krg.getProcessNode(head.process.type)
              return [
                {
                  id: head.process.id,
                  kind: 'process' as 'process',
                  label: process.meta.label,
                  color: head.id === node_id ? '#B3CFFF' : 'lightgrey',
                  icon: process.meta.icon || [func_icon],
                  parents: Object.keys(head.process.inputs).length === 0 ? ['start'] : Object.values(head.process.inputs).map(({ id }) => `${id}:output`),
                },
                {
                  id: `${head.process.id}:output`,
                  kind: 'data' as 'data',
                  label: process.output.meta.label,
                  color: head.id === node_id ? '#B3CFFF' : 'lightgrey',
                  icon: process.output.meta.icon || [variable_icon],
                  parents: [head.process.id],
                },
              ]
            }),
            {
              id: 'extend',
              kind: 'process' as 'process',
              label: 'Extend',
              color: extend || suggest ? '#B3CFFF' : 'lightgrey',
              content: '+',
              parents: [head ? `${head.process.id}:output` : 'start'],
            },
          ]}
          onclick={(_evt, id) => {
            if (id === 'extend') {
              router.push(`/graph/${graph_id}${graph_id !== node_id ? `/node/${node_id}` : ''}/extend`, undefined, { shallow: true })
            } else {
              const focus_node_process_id = id.split(':')[0]
              const focus_node_id = id === 'start' ? 'start' : metapath.filter((node) => node.process.id === focus_node_process_id)[0].id
              router.push(`/graph/${graph_id}${graph_id !== focus_node_id ? `/node/${focus_node_id}` : ''}`, undefined, { shallow: true })
            }
          }}
        />
      </div>
      <main className="flex-grow flex flex-col">
        {error ? <div>{error}</div> : null}
        {suggest ?
          <Suggest id={graph_id} head={head} />
          : extend ?
            <Extend id={graph_id} head={head} />
            : node_id === 'start' ?
              <Home />
              : head ?
                <Cell id={graph_id} head={head} />
                : null}
      </main>
    </>
  )
}

