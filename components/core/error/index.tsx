import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { error_icon } from '@/icons'

export const Error = MetaNode('Error')
  .meta({
    label: 'Error',
    description: 'An error has occurred',
    icon: [error_icon],
    example: 'An unknown error occurred',
  })
  .codec(z.string())
  .view(error => (
    <div>{error}</div>
  ))
  .build()
