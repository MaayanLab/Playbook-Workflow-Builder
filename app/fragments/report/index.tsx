import React from 'react'
import { useExRouter } from '@/app/fragments/ex-router'
import { z } from 'zod'
import dynamic from 'next/dynamic'
import { MetapathProvider, useKRG } from '@/app/fragments/metapath'
import { usePathname } from 'next/navigation'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const Cells = dynamic(() => import('@/app/fragments/report/cells'))
const Chat = dynamic(() => import('@/app/fragments/chat/chat'))

const QueryType = z.object({
  id: z.string().optional(),
  session_id: z.string().optional(),
  thread_id: z.string().optional(),
})

function Report({ params }: { params: z.infer<typeof QueryType> }) {
  const krg = useKRG()
  return (
    <Layout
      sidebar={<Chat mode="report" graph_id={params.id} embedded />}
    >
      <main className="grow flex flex-col mx-4">
        {params.id ?
          <Cells
            session_id={params.session_id}
            krg={krg}
            id={params.id}
          />
          : <div className="alert alert-error">Page not found</div>}
      </main>
    </Layout>
  )
}

export default function App({ routeParams }: { routeParams?: z.infer<typeof QueryType> } = {}) {
  const router = useExRouter()
  const pathname = usePathname()
  // Input edits update browser history without remounting this report page.
  // Route props retain their server value until the next server navigation.
  const params = React.useMemo(() => {
    const match = pathname?.match(/^(?:\/session\/([^/]+))?\/report(?:\/([^/]+))?\/?$/)
    return QueryType.parse(match
      ? { ...(routeParams ?? router.query), session_id: match[1], id: match[2] }
      : routeParams ?? router.query)
  }, [pathname, routeParams, router.query])
  return (
    <MetapathProvider session_id={params.session_id}>
      <Report params={params} />
    </MetapathProvider>
  )
}
