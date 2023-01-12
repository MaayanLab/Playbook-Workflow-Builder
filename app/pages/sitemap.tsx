import React from 'react'
import Head from 'next/head'
import dynamic from 'next/dynamic'
import { Card, Elevation } from "@blueprintjs/core";
import Link from 'next/link';
 
const Header = dynamic(() => import('@/app/fragments/playbook/header'))
const Footer = dynamic(() => import('@/app/fragments/playbook/footer'))
const Contributors = dynamic(() => import('@/app/fragments/playbook/contributors'))

export default function SiteMap() {
  return (
    <div className="flex flex-col min-w-screen min-h-screen">
      <Head>
        <title>Playbook</title>
      </Head>
      <Header homepage="/" />
        <main className="flex-grow container mx-auto py-4 flex flex-col">
          <h2 className="bp4-heading">Site Map</h2>
          <p className="bp4-text-large">Click the card to go to the relevant page.</p>
          <div className="grid md:grid-cols-3 gap-4 my-2">
            <Link href="/graph">
              <Card className="col-span-1" interactive={true} elevation={Elevation.TWO}>
                <h5 className="bp4-heading">Graph Mode</h5>
                <p>A graph interface for incrementally building complex workflows of MetaNodes.</p>
              </Card>
            </Link>
            <Link href="/report">
              <Card className="col-span-1" interactive={true} elevation={Elevation.TWO}>
                <h5 className="bp4-heading">Report Mode</h5>
                <p>A jupyter-notebook style interface for sequential execution of MetaNode Processes.</p>
              </Card>
            </Link>
            <Link href="/testing">
              <Card className="col-span-1" interactive={true} elevation={Elevation.TWO}>
                <h5 className="bp4-heading">Prototyping Interface</h5>
                <p>A simple interface for developing new components & inspecting internals.</p>
              </Card>
            </Link>
          </div>
          <Contributors />
        </main>
      <Footer />
    </div>
  )
}