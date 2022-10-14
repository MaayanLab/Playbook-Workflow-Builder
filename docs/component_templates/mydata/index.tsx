import React from 'react'
import { MetaNode } from '@/spec/metanode'
import * as t from 'io-ts'
import codecFrom from '@/utils/io-ts-codec'

// A unique name for your data type is used here
export const MyData = MetaNode.createData('MyData')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'My Data',
    description: 'My data type, rendered',
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using io-ts, a compile-time and runtime type-safe codec can be constructed
  .codec(codecFrom(t.type({
    field: t.array(
      t.partial({
        maybe: t.string
      })
    ),
  })))
  // react component rendering your data goes here
  .view(data => {
    return (
      <div>{JSON.stringify(data)}</div>
    )
  })
  .build()
