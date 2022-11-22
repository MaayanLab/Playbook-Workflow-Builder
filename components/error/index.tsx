import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import codecFrom from '@/utils/zod-codec'

export const Error = MetaNode.createData('Error')
  .meta({
    label: 'Error',
    description: 'An error has occurred',
  })
  .codec(codecFrom(z.string()))
  .view(error => (
    <div>{error}</div>
  ))
  .build()
