import Link from 'next/link'
import dynamic from 'next/dynamic'
import krg from "@/app/krg"
import { func_icon, variable_icon } from '@/icons'
import { MetaNode, DataMetaNode, ProcessMetaNode } from '@/spec/metanode'

const Markdown = dynamic(() => import('@/app/components/Markdown'))
const Icon = dynamic(() => import('@/app/components/icon'))

export function metanodeType(metanode: MetaNode) {
  if (metanode.kind === 'data') {
    const terminal = krg.getNextProcess(metanode.spec).length === 0
    if (terminal) return 'Viewer'
    else return 'Data'
  } else {
    const terminal = krg.getNextProcess(metanode.output.spec).length === 0
    if ('prompt' in metanode) {
      if (terminal) return 'Interactive Viewer'
      else return 'Prompt'
    } else {
      if (terminal) return 'View Resolver'
      else return 'Resolver'
    }
  }
}

export default function MetaNodeListing({ metanodes }: { metanodes: { metanode: (DataMetaNode | ProcessMetaNode), multi?: boolean, type?: string, }[] }) {
  return (
    <div className="overflow-x-auto">
      <table className="table">
        <thead>
          <tr>
            <th className="prose max-w-none">Icon</th>
            <th className="prose max-w-none">Type</th>
            <th className="prose max-w-none">Label</th>
            <th className="prose max-w-none">Description</th>
          </tr>
        </thead>
        <tbody>
          {metanodes.map(({ metanode, type }) => {
            return (
              <tr key={`${type}-${metanode.spec}`}>
                <td><Icon icon={metanode.meta.icon ? metanode.meta.icon : metanode.kind === 'data' ? variable_icon : func_icon} /></td>
                <td>{type ? type : metanodeType(metanode)}</td>
                <th>
                  <Link
                    href={`/components/${encodeURIComponent(metanode.spec)}`}
                    className="link prose max-w-none whitespace-pre-wrap"
                  >{metanode.meta.label}</Link>
                </th>
                <td className="prose max-w-none whitespace-pre-wrap">
                  <Markdown components={{
                    p(props: React.PropsWithChildren<{}>) {
                      return <span {...props}>{props.children}</span>
                    }
                  }}>{metanode.meta.description}</Markdown>
                </td>
              </tr>
            )
          })}
        </tbody>
      </table>
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
