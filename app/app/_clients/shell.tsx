'use client'

import React from 'react'
import dynamic from 'next/dynamic'
import { usePathname, useSearchParams } from 'next/navigation'
import { MetapathProvider } from '@/app/fragments/metapath'
import { useIsClient } from '@/app/fragments/client-only'
import { parsePathname } from '@/app/app/_lib/routes'

// only rendered in the browser
const SessionProvider = dynamic(() => import('@/app/extensions/auth-client').then(module => module.SessionProvider))
const HotkeysProvider = dynamic(() => import('@blueprintjs/core').then(module => module.HotkeysProvider))
const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const Chat = dynamic(() => import('@/app/fragments/chat/chat'))

/**
 * See @/app/app/(shell)/layout: the playbook state (MetapathProvider, e.g. what's known about each
 *  process' resolution) is shared by the pages, and the chat sidebar of the graph & report pages
 *  is the same chat while moving between them.
 *
 * On the server (& while hydrating) only the page is rendered (its browser only views render
 *  nothing there, see ClientOnly): rendering it, rather than nothing, has the response wait for the
 *  page, so that its `notFound()` still responds with a 404 status.
 */
export default function Shell({ children }: React.PropsWithChildren) {
  const isClient = useIsClient()
  const pathname = usePathname() ?? ''
  const searchParams = useSearchParams()
  const route = React.useMemo(() => parsePathname(pathname), [pathname])
  const thread_id = searchParams?.get('thread_id') ?? undefined
  if (!isClient) return <>{children}</>
  const sidebar = route.kind === 'graph' ? <Chat mode="graph" graph_id={route.graph_id} node_id={route.node_id ?? route.graph_id} embedded />
    : route.kind === 'report' ? <Chat mode="report" graph_id={route.id} thread_id={thread_id} embedded />
    : undefined
  return (
    <SessionProvider>
      <HotkeysProvider>
        <MetapathProvider session_id={route.session_id}>
          <Layout sidebar={sidebar}>{children}</Layout>
        </MetapathProvider>
      </HotkeysProvider>
    </SessionProvider>
  )
}
