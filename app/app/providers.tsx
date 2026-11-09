'use client'

import dynamic from 'next/dynamic'
import type { PropsWithChildren } from 'react'
import { TRPCReactProvider } from '@/app/trpc/client'
import { RuntimeConfig } from '@/app/fragments/config'

const Analytics = dynamic(() => import('@/app/fragments/analytics'), { ssr: false })

/**
 * Providers which can render on the server, browser-only parts are in
 *  @/app/app/_clients/client-page (rendered by each page)
 */
export default function Providers({ runtimeConfig, children }: PropsWithChildren<{ runtimeConfig: Record<string, string> }>) {
  return (
    <TRPCReactProvider>
      <RuntimeConfig value={runtimeConfig}>
        {children}
        <Analytics />
      </RuntimeConfig>
    </TRPCReactProvider>
  )
}
