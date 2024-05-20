import '@blueprintjs/icons/lib/css/blueprint-icons.css'
import '@blueprintjs/select/lib/css/blueprint-select.css'
import '@blueprintjs/core/lib/css/blueprint.css'
import '@blueprintjs/table/lib/css/table.css'
import '@/app/styles/styles.css'
import type { AppProps } from 'next/app'
import Head from 'next/head'
import dynamic from 'next/dynamic'

const SessionProvider = dynamic(() => import('next-auth/react').then(({ SessionProvider }) => SessionProvider))
const HotkeysProvider = dynamic(() => import('@blueprintjs/core').then(({ HotkeysProvider }) => HotkeysProvider))
const RuntimeConfig = dynamic(() => import('@/app/fragments/config').then(({ RuntimeConfig }) => RuntimeConfig), { ssr: false })
const Analytics = dynamic(() => import('@/app/fragments/analytics'), { ssr: false })

export default function App({ Component, pageProps: { session, ...pageProps } }: AppProps & { pageProps: { session: any } }) {
  return (
    <>
      <Head>
        <title>Playbook</title>
        <meta name="viewport" content="width=device-width, initial-scale=1.0" />
        <link rel="icon" type="image/x-icon" href="/favicon.ico" />
        <link rel="shortcut icon" type="image/x-icon" href="/favicon.ico" />
      </Head>
      <RuntimeConfig>
        <SessionProvider session={session}>
          <HotkeysProvider>
            <Component {...pageProps} />
          </HotkeysProvider>
        </SessionProvider>
        <Analytics />
      </RuntimeConfig>
    </>
  )
}
