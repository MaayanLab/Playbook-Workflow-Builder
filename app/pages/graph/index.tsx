import '@blueprintjs/icons/lib/css/blueprint-icons.css'
import '@blueprintjs/select/lib/css/blueprint-select.css'
import '@blueprintjs/core/lib/css/blueprint.css'

import dynamic from 'next/dynamic'
import Head from 'next/head'
import { useRouter } from 'next/router'

import { start_icon } from '@/icons'

const Header = dynamic(() => import('@/app/fragments/playbook/header'))
const Home = dynamic(() => import('@/app/fragments/playbook/home'))
const Footer = dynamic(() => import('@/app/fragments/playbook/footer'))
const Breadcrumbs = dynamic(() => import('@/app/fragments/breadcrumbs'))

export default function App() {
  const router = useRouter()
  return (
    <div className="flex flex-col min-w-screen min-h-screen">
      <Head>
        <title>Playbook</title>
      </Head>

      <Header />

      <div className="flex w-auto h-40">
        <Breadcrumbs
          graph={[
            {
              id: 'start',
              kind: 'data',
              label: 'Start',
              color: '#B3CFFF',
              icon: [start_icon],
              parents: [],
            },
            {
              id: 'extend',
              kind: 'process',
              label: 'Extend',
              color: 'lightgrey',
              content: '+',
              parents: ['start'],
            }
          ]}
          onclick={(_evt, id) => {
            if (id === 'extend') {
              router.push(`/graph/extend`, undefined, { shallow: true })
            }
          }}
        />
      </div>

      <main className="flex-grow m-4">
        <Home />
      </main>

      <Footer />
    </div>
  )
}