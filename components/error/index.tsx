import React from 'react'
import { MetaNode } from '@/spec/metanode'
import * as t from 'io-ts'
import codecFrom from '@/utils/io-ts-codec'

export const Error = MetaNode.createData('Error')
  .meta({
    label: 'Error',
    description: 'An error has occurred',
  })
  .codec(codecFrom(t.string))
  .view(error => (
    <div>{error}</div>
  ))
  .build()
