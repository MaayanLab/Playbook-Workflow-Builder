import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'

export const Error = MetaNode.createData('Error')
  .meta({
    label: 'Error',
    description: 'An error has occurred',
  })
  .codec(z.string())
  .view(error => (
    <div>{error}</div>
  ))
  .build()
