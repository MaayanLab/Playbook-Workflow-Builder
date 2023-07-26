import React from 'react'
import { z } from 'zod'
import { MetaNode } from "@/spec/metanode"
import { AnnData } from '@/components/data/anndata'
import { useClientMetadataFromFile } from './api/metadata/client'
import { useClientUpdateMetadata } from './api/metadata/update/client'

export const Label = MetaNode('Label')
  .meta({
    label: 'Label',
    description:' label',
  })
  .inputs({ matrix: AnnData })
  .output(AnnData)
  .prompt(props => {
    const { data } = useClientMetadataFromFile(props.inputs.matrix)
    const { trigger } = useClientUpdateMetadata(props.inputs.matrix.url)
    return (
      <>
      {JSON.stringify(data)}
      <button onClick={async () => {
        // here I just copy condition into condition condition2 column
        //  --instead of this you want to specify some known column that will
        //  be used by differential expression modules
        props.submit(
          await trigger({
            file: props.inputs.matrix,
            data: {
              'condition2': (data as any)['condition'],
            },
          })
        )
      }}>Submit</button>
      </>
    )
  })
  .story(props => ``)
  .build()
