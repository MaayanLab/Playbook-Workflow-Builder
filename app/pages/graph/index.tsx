import '@blueprintjs/icons/lib/css/blueprint-icons.css'
import '@blueprintjs/select/lib/css/blueprint-select.css'
import '@blueprintjs/core/lib/css/blueprint.css'

import dynamic from 'next/dynamic'
import Head from 'next/head'
import { useRouter } from 'next/router'

import { start_icon } from '@/icons'

// const Viz = dynamic(() => import('@/app/fragments/metagraph-viz'))
const Viz = () => null
const Header = dynamic(() => import('@/app/fragments/header'))
const Footer = dynamic(() => import('@/app/fragments/footer'))
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
        <div className="bp4-running-text">
          <h1 className="bp4-heading">Welcome to Signature Commons Graph</h1>
          <p className="bp4-text-large bp4-text-muted">An interactive knowledge enrichment tool for user-driven meta-graph traversal</p>
          <hr />
          <p>
            You can start <i>extending</i> your own instance graph by clicking the <svg className="inline-block" width={32} height={32} viewBox="0 0 1 1"><rect width={1} height={1} fill="lightgrey"></rect><text x={0.5} y={0.5} fontSize="0.5px" dominantBaseline="middle" textAnchor="middle" fill="black">+</text></svg> above.
          </p>
          <p>
            To completely start over rather than expanding from Home, click the Signature Commons Graph logo.
          </p>
          <figure>
            <Viz />
            <figcaption className="text-center">A global view of the meta graph</figcaption>
          </figure>
        </div>
      </main>

      <Footer />
    </div>
  )
}