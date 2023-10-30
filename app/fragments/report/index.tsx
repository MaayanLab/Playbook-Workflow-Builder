import React from 'react'
import type { GetServerSidePropsContext } from 'next'
import { useRouter } from 'next/router'
import fpprg from '@/app/fpprg'
import krg from '@/app/krg'
import { z } from 'zod'
import { SWRConfig } from 'swr'
import dynamic from 'next/dynamic'
import fetcher from '@/utils/next-rest-fetcher'
import { MetapathProvider } from '@/app/fragments/metapath'

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const Cells = dynamic(() => import('@/app/fragments/report/cells'))

const QueryType = z.object({
  id: z.string().optional(),
  session_id: z.string().optional(),
})

export async function getServerSideProps(ctx: GetServerSidePropsContext) {
  const params = QueryType.parse(ctx.params || {})
  if (!params.id || params.session_id !== undefined) {
    return {
      props: {
        fallback: {}
      }
    }
  }
  const fpl = await fpprg.getFPL(params.id)
  if (fpl === undefined) {
    return { notFound: true }
  }
  const results = fpl.resolve().map(fpl => fpl.toJSON())
  const fallback: Record<string, unknown> = {
    [`/api/db/fpl/${params.id}`]: results
  }
  for (const result of results) {
    const output = await fpprg.getResolved(result.process.id)
    if (output) fallback[`/api/db/process/${result.process.id}/output`] = output.toJSON().data
  }
  return {
    props: {
      fallback: JSON.parse(JSON.stringify(fallback)),
    }
  }
}

export default function App({ fallback }: { fallback: any }) {
  const router = useRouter()
  const params = QueryType.parse(router.query)
  return (
    <Layout>
      <SWRConfig value={{ fallback, fetcher }}>
        <main className="flex-grow container mx-auto py-4 flex flex-col">
          {params.id ? 
            <MetapathProvider
              session_id={params?.session_id}
              id={params.id}
            >
              <Cells
                session_id={params?.session_id}
                krg={krg}
                id={params.id}
              />
            </MetapathProvider>
          : <div className="alert alert-error">Page not found</div>}
        </main>
      </SWRConfig>
    </Layout>
  )
}
