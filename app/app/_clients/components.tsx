'use client'

import React from 'react'
import krg from "@/app/krg"
import dynamic from 'next/dynamic'

const MetaNodeListing = dynamic(() => import('@/app/fragments/components/listing'))

export default function ComponentPage() {
  return (
    <>
      <div className="container flex-grow self-center flex place-content-center overflow-hidden py-5">
        <MetaNodeListing metanodes={[...krg.getDataNodes(), ...krg.getProcessNodes()].map(metanode => ({ metanode }))} />
      </div>
    </>
  )
}
