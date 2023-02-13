import '@blueprintjs/icons/lib/css/blueprint-icons.css'
import '@blueprintjs/core/lib/css/blueprint.css'
import '@blueprintjs/table/lib/css/table.css'
import '@blueprintjs/popover2/lib/css/blueprint-popover2.css'
import '@/app/styles/styles.css'
import type { AppProps } from 'next/app'
import { HotkeysProvider } from '@blueprintjs/core'
import { SessionProvider } from 'next-auth/react'
import dynamic from 'next/dynamic'

const Analytics = dynamic(() => import('@/app/fragments/analytics'), { ssr: false })

export default function App({ Component, pageProps: { session, ...pageProps } }: AppProps & { pageProps: { session: any } }) {
  return (
    <>
      <SessionProvider session={session}><HotkeysProvider><Component {...pageProps} /></HotkeysProvider></SessionProvider>
      <Analytics />
    </>
  )
}
