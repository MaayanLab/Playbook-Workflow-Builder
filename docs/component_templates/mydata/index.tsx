import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'

// A unique name for your data type is used here
export const MyData = MetaNode('MyData')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'My Data',
    description: 'My data type, rendered',
  })
  // this should have a codec which can encode or decode the data type represented by this node
  //  using zod, a compile-time and runtime type-safe codec can be constructed
  .codec(z.object({
    field: z.array(
      z.object({
        maybe: z.string().optional()
      })
    ),
  }))
  // react component rendering your data goes here
  .view(data => {
    return (
      <div>{JSON.stringify(data)}</div>
    )
  })
  .build()
