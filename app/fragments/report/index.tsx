import React from 'react'
import type { GetServerSidePropsContext } from 'next'
import { useExRouter } from '@/app/fragments/ex-router'
import fpprg from '@/app/fpprg'
import useKRG from '@/app/fragments/graph/krg'
import { z } from 'zod'
import dynamic from 'next/dynamic'
import { MetapathProvider } from '@/app/fragments/metapath'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const Cells = dynamic(() => import('@/app/fragments/report/cells'))
const Chat = dynamic(() => import('@/app/fragments/chat/chat'))

const QueryType = z.object({
  id: z.string().optional(),
  session_id: z.string().optional(),
  thread_id: z.string().optional(),
})

export async function getServerSideProps(ctx: GetServerSidePropsContext) {
  const params = QueryType.parse(ctx.params || {})
  if (!params.id || params.session_id !== undefined) {
    return {
      props: {}
    }
  }
  const fpl = await fpprg.getFPL(params.id)
  if (fpl === undefined) {
    return { notFound: true }
  }
  const results = fpl.resolve().map(fpl => fpl.toJSON())
  return { props: {} }
}

export default function App() {
  const router = useExRouter()
  const params = React.useMemo(() => QueryType.parse(router.query), [router.query])
  const krg = useKRG({ session_id: params?.session_id })
  return (
    <>
      <MetapathProvider session_id={params?.session_id}>
        <Layout
          sidebar={<Chat mode="report" graph_id={params.id} embedded />}
        >
          <main className="grow flex flex-col mx-4">
            {params.id ? 
              <Cells
                session_id={params?.session_id}
                krg={krg}
                id={params.id}
              />
              : <div className="alert alert-error">Page not found</div>}
          </main>
        </Layout>
      </MetapathProvider>
    </>
  )
}
