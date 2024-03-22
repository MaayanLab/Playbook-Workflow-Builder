import React from 'react'
import Head from 'next/head'
import dynamic from 'next/dynamic'
import Link from 'next/link';

const Layout = dynamic(() => import('@/app/fragments/playbook/layout'))
const Contributors = dynamic(() => import('@/app/fragments/playbook/contributors'))

export default function SiteMap() {
  return (
    <Layout>
      <Head>
        <title>Playbook: Site Map</title>
      </Head>
      <main className="flex-grow container mx-auto py-4 flex flex-col">
        <div className="hero">
          <div className="hero-content prose max-w-none">
            <div>
              <h2>Site Map</h2>
              <p>Click the card to go to the relevant page.</p>
              <div className="grid md:grid-cols-3 gap-4 my-2">
                <Link href="/graph/extend">
                  <div className="card card-compact card-bordered col-span-1 bg-base-100 shadow shadow-slate-300 hover:shadow-slate-600 hover:cursor-pointer">
                    <div className="card-body prose max-w-none">
                      <h4 className="card-title">Graph Mode</h4>
                      <p>A graph interface for incrementally building complex workflows of MetaNodes.</p>
                    </div>
                  </div>
                </Link>
                <Link href="/playbooks">
                  <div className="card card-compact card-bordered col-span-1 bg-base-100 shadow shadow-slate-300 hover:shadow-slate-600 hover:cursor-pointer">
                    <div className="card-body prose max-w-none">
                      <h4 className="card-title">Playbook Catalog</h4>
                      <p>An interface for searching for and launching published playbooks.</p>
                    </div>
                  </div>
                </Link>
                <Link href="/testing">
                  <div className="card card-compact card-bordered col-span-1 bg-base-100 shadow shadow-slate-300 hover:shadow-slate-600 hover:cursor-pointer">
                    <div className="card-body prose max-w-none">
                      <h4 className="card-title">Prototyping Interface</h4>
                      <p>A simple interface for developing new components & inspecting internals.</p>
                    </div>
                  </div>
                </Link>
                <Link href="/explore">
                  <div className="card card-compact card-bordered col-span-1 bg-base-100 shadow shadow-slate-300 hover:shadow-slate-600 hover:cursor-pointer">
                    <div className="card-body prose max-w-none">
                      <h4 className="card-title">Exploration Visualization</h4>
                      <p>An interactive d3 playground for exploring the metagraph.</p>
                    </div>
                  </div>
                </Link>
              </div>
            </div>
          </div>
        </div>
        <Contributors />
      </main>
    </Layout>
  )
}