import '@blueprintjs/icons/lib/css/blueprint-icons.css'
import '@blueprintjs/select/lib/css/blueprint-select.css'
import '@blueprintjs/core/lib/css/blueprint.css'
import '@blueprintjs/table/lib/css/table.css'
import '@/app/styles/styles.css'
import type { Metadata } from 'next'
import type { PropsWithChildren } from 'react'
import Providers from './providers'
import { getRuntimeConfig } from '@/app/runtime-config'

// the runtime configuration is read from the environment on each request (not at build time)
export const dynamic = 'force-dynamic'

export const metadata: Metadata = {
  title: 'Playbook',
  icons: { icon: '/favicon.ico', shortcut: '/favicon.ico' },
}

export default function RootLayout({ children }: PropsWithChildren) {
  return (
    <html lang="en">
      <body><Providers runtimeConfig={getRuntimeConfig()}>{children}</Providers></body>
    </html>
  )
}
