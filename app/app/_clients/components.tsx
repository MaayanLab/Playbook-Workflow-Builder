'use client'

import React from 'react'
import krg from "@/app/krg"
import dynamic from 'next/dynamic'

const MetaNodeListing = dynamic(() => import('@/app/fragments/components/listing'))

export default function ComponentPage() {
  return (
    <>
      {/* (no overflow clipping here: the listing's toolbar sticks to the top of the page's scroll) */}
      <div className="container mx-auto px-4 py-5 flex flex-col gap-2">
        <h1 className="prose max-w-none text-3xl font-bold">Component Catalog</h1>
        <MetaNodeListing metanodes={[...krg.getDataNodes(), ...krg.getProcessNodes()].map(metanode => ({ metanode }))} />
      </div>
    </>
  )
}
