'use client'

import dynamic from 'next/dynamic'
import type { PropsWithChildren } from 'react'
import ClientOnly from '@/app/fragments/client-only'

// only rendered in the browser (inside ClientOnly)
const SessionProvider = dynamic(() => import('@/app/extensions/auth-client').then(module => module.SessionProvider))
const HotkeysProvider = dynamic(() => import('@blueprintjs/core').then(module => module.HotkeysProvider))

/**
 * The client side of a page: rendered in the browser only. Pages render this
 *  rather than the root layout so that the page itself still renders on the server
 *  (e.g. `notFound()` responds with a 404 status).
 */
export default function ClientPage({ children }: PropsWithChildren) {
  return (
    <ClientOnly>
      <SessionProvider>
        <HotkeysProvider>{children}</HotkeysProvider>
      </SessionProvider>
    </ClientOnly>
  )
}
