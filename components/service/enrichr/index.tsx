import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { GeneSet } from '@/components/core/input/set'
import { z } from 'zod'
import { gene_icon, enrichr_icon } from '@/icons'

const enrichr_url = 'https://maayanlab.cloud/Enrichr'

export const EnrichrUserList = MetaNode.createData('EnrichrUserList')
  .meta({
    label: 'Enrichr User List',
    description: 'A gene set submitted to Enrichr',
    icon: [enrichr_icon, gene_icon],
  })
  .codec(z.object({
    shortId: z.string(),
    userListId: z.number(),
  }))
  .view(userlist => (
    <div className="flex-grow flex flex-row m-0" style={{ minHeight: 500 }}>
      <iframe
        className="flex-grow border-0"
        src={`${enrichr_url}/enrich?dataset=${userlist.shortId}`}
      />
    </div>
  ))
  .build()

export const EnrichrGenesetSearch = MetaNode.createProcess('EnrichrGenesetSearch')
  .meta({
    label: `Enrichr Enrichment Analysis`,
    tags: {
      'Input Type': {
        Gene: 1,
      },
      'Input Cardinality': {
        Set: 1,
      },
      Service: {
        Enrichr: 1,
      },
      'Output Type': {
        Signature: 1,
        Gene: 1,
      },
      'Output Cardinality': {
        Weighted: 1,
      },
    },
    icon: [enrichr_icon],
    description: "Fisher's exact test, odd ratio, Jaccard index",
  })
  .inputs({ geneset: GeneSet })
  .output(EnrichrUserList)
  .resolve(async (props) => {
    const formData = new FormData()
    formData.append('list', props.inputs.geneset.join('\n'))
    formData.append('description', `playbook-partnership`)
    const response = await fetch(`${enrichr_url}/addList`, {
      method: 'post',
      body: formData,
    })
    if (response.status !== 200) {
      throw new Error(`'Enrichr addList' Request failed with status ${response.status}: ${await response.text()}`)
    }
    return await response.json()
  })
  .build()
