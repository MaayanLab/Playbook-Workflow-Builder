import React from 'react'
import krg from "@/app/krg"
import { useRouter } from 'next/router'
import dynamic from 'next/dynamic'
import Link from 'next/link'
import * as dict from '@/utils/dict'
import { func_icon, variable_icon } from '@/icons'
import MetaNodeListing, { metanodeType } from '@/app/fragments/components/listing'

const Markdown = dynamic(() => import('@/app/components/Markdown'))
const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const Icon = dynamic(() => import('@/app/components/icon'))

export function getServerSideProps({ query }: { query: { spec: string } }) {
  const metanode = typeof query.spec === 'string' ? krg.getDataNode(query.spec) ?? krg.getProcessNode(query.spec) : undefined
  if (metanode === undefined) return { notFound: true }
  return {
    props: {
      spec: metanode.spec,
      kind: metanode.kind,
    },
  }
}

export default function ComponentPage({ spec, kind }: { spec: string, kind: 'data' | 'process' }) {
  const router = useRouter()
  const { metanode, prevMetanodes, nextMetanodes, inputMetanodes, outputMetanode } = React.useMemo(() => {
    const dataMetanode = typeof router.query.spec === 'string' ? krg.getDataNode(spec) : undefined
    const processMetanode = typeof router.query.spec === 'string' ? krg.getProcessNode(spec) : undefined
    const metanode = dataMetanode ? dataMetanode : processMetanode
    const nextMetanodes = dataMetanode ? krg.getNextProcess(dataMetanode.spec) : undefined
    const prevMetanodes = dataMetanode ? krg.getPrevProcess(dataMetanode.spec) : undefined
    const inputMetanodes = processMetanode && processMetanode.inputs !== undefined && !dict.isEmpty(processMetanode.inputs) ? processMetanode.inputs : undefined
    const outputMetanode = processMetanode ? processMetanode.output : undefined
    return { metanode, prevMetanodes, nextMetanodes, inputMetanodes, outputMetanode }
  }, [router.query.spec])
  if (metanode === undefined) return <div>Not found</div>
  return (
    <Layout>
      <div className="container flex-grow self-center flex flex-col py-5 overflow-hidden">
        <div className="flex flex-row items-center gap-2 m-2">
          <Icon icon={metanode.meta.icon ? metanode.meta.icon : metanode.kind === 'data' ? variable_icon : func_icon} size={3} title={null} />
          <div className="flex flex-col prose max-w-none">
            <h1 className="m-0 whitespace-pre-wrap">{metanode.meta.label}</h1>
            <h2 className="m-0 text-gray-500">{metanodeType(metanode)}</h2>
          </div>
        </div>
        {metanode.meta.version ? <p className="prose max-w-none"><strong>Version</strong>: {metanode.meta.version}</p> : null}
        {metanode.meta.author ? <p className="prose max-w-none"><strong>Author</strong>: {metanode.meta.author}</p> : null}
        {metanode.meta.license ? <p className="prose max-w-none"><strong>License</strong>: {metanode.meta.license}</p> : null}
        {metanode.meta.tags ? <p className="prose max-w-none"><strong>Tags</strong>: <div className="inline-flex flex-row gap-2 flex-wrap">{dict.items(metanode.meta.tags).flatMap(({ key, value }) => dict.keys(value).map(v => `${v} (${key})`)).map(tag =>
          <div className="badge badge-secondary prose max-w-none" key={tag}>{tag}</div>
        )}</div>
        </p> : null}
        <Markdown components={{
          p(props: React.PropsWithChildren<{}>) {
            return <p className="prose max-w-none" {...props}>{props.children}</p>
          }
        }}>{`**Description**: ${metanode.meta.description}`}</Markdown>
          <MetaNodeListing
            metanodes={[
              ...inputMetanodes ? dict.items(inputMetanodes)
                .flatMap(({ key, value }) => (Array.isArray(value) ? value : [value]).map(v =>
                  ({ metanode: v, type: Array.isArray(value) ? `Multi-Input Data` : 'Input Data' })
                )) : [],
              ...outputMetanode ? [{ metanode: outputMetanode, type: `Output ${metanodeType(outputMetanode)}` }] : [],
              ...prevMetanodes ? prevMetanodes.map(metanode => ({ metanode, type: `Previous ${metanodeType(metanode)}` })) : [],
              ...nextMetanodes ? nextMetanodes.map(metanode => ({ metanode, type: `Next ${metanodeType(metanode)}` })) : [],
            ]}
          />
        <div className="justify-self-end"><Link href="/components" className="btn btn-sm">&lt;- All Components</Link></div>
      </div>
    </Layout>
  )
}
