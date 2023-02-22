import krg from '@/app/krg'
import dynamic from 'next/dynamic'
import { start_icon, variable_icon } from '@/icons'
import Link from 'next/link'

const Icon = dynamic(() => import('@/app/components/icon'))
const Card = dynamic(() => import('@blueprintjs/core').then(({ Card }) => Card))

export default function ChooseStart() {
  return (
    <div className="container mx-auto flex flex-col items-center py-4">
      <h1 className="bp4-heading"><Icon icon={start_icon} />&nbsp;Start</h1>
      <div className="flex flex-row flex-wrap gap-2 justify-center">
        {krg.getDataNodes().map(item => {
          return (
            <Link href={`/explore/dfs/${item.spec}`}>
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
        })}
      </div>
    </div>
  )
}