import React from 'react'
import { MetaNode } from '@/spec/metanode'

export type GlyGenDataType = {
  mydatatype: string[]
}

export type GlyGenResponseType = {
  mydatatype: string[]
}

export const GlyGenResponse = MetaNode.createData('GlyGenResponse')
  .meta({
    label: 'GlyGen Response',
    description: 'GlyGen response object',
  })
  .codec<GlyGenResponseType>()
  .view(data => (
    <div>{JSON.stringify(data)}</div>
  ))
  .build()

export const GlyGenQuery = MetaNode.createData('GlyGenQuery')
  .meta({
    label: 'GlyGen Query',
    description: 'Search for records in GlyGen. Must be a JSON input',
  })
  .codec<GlyGenDataType>()
  .view(data => (
    <div>{JSON.stringify(data)}</div>
  ))
  .build()

export const ProteinProductInformation = MetaNode.createProcess('ProteinProductInformation')
  .meta({
    label: 'Protein Product Information',
    description: 'Search for protein records in GlyGen',
  })
  .inputs({ query: GlyGenQuery })
  .output(GlyGenResponse)
  .resolve(async (props) => {
    const query = JSON.stringify(props.inputs.query)
    const response = await fetch(`https://api.glygen.org/directsearch/protein/?query=${query}`, {
      method: 'GET',
      headers: {
        accept: 'application/json'
      }
      })
    return response.json()
  })
  .build()