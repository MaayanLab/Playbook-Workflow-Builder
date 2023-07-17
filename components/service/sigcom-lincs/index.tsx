import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { ScoredDrugs, ScoredGenes } from '@/components/core/input/scored'

const sigcom_lincs_url = 'https://maayanlab.cloud/sigcom-lincs'

export const SigComLINCSSignatureResults = MetaNode(`SigComLINCSSignatureResults`)
  .meta({
    label: 'SigCom LINCS Signature Search Results',
    description: 'LINCS L1000 Signatures Query Results',
  })
  .codec(z.object({ id: z.string() }))
  .view(props => (
    <div className="flex-grow flex flex-row m-0" style={{ minHeight: 800 }}>
      <iframe
        className="flex-grow border-0"
        src={`${sigcom_lincs_url}/#/SignatureSearch/Rank/${props.id}`}
      />
    </div>
  ))
  .build()

export const SigComLINCSSignatureSearch = MetaNode(`SigComLINCSSignatureSearch`)
  .meta({
    label: 'SigCom LINCS Signature Search',
    description: 'Query LINCS L1000 Signatures',
  })
  .inputs({ genes: ScoredGenes })
  .output(SigComLINCSSignatureResults)
  .resolve(async (props) => {
    const up: Record<string, true> = {}
    const down: Record<string, true> = {}
    props.inputs.genes.forEach(({ term, zscore }) => {
      if (typeof zscore !== 'number') return
      if (zscore >= 3) up[term] = true
      else if (zscore <= 3) down[term] = true
    })
    const metaEntityReq = await fetch(`${sigcom_lincs_url}/metadata-api/entities/find`, {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        filter: {
          where: {
            'meta.symbol': {
              inq: [...Object.keys(up), ...Object.keys(down)],
            }
          }
        }
      }),
    })
    const metaRes = await metaEntityReq.json()
    const up_entities: string[] = []
    const down_entities: string[] = []
    metaRes.forEach((ent: any) => {
      if (up[ent.meta.symbol]) up_entities.push(ent.id)
      if (down[ent.meta.symbol]) down_entities.push(ent.id)
    })
    const metaInputReq = await fetch(`${sigcom_lincs_url}/metadata-api/user_input`, {
      method: 'POST',
      headers: {
        'Accept': 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        meta: {
            "$validator": "/dcic/signature-commons-schema/v6/meta/user_input/user_input.json",
            up_entities,
            down_entities,
        },
        type: "signature",
      }),
    })
    const results = await metaInputReq.json()
    return { id: results.id }
  })
  .story(props => `Identified reversers and mimickers from over 1 million signatures using SigCom LINCS [\\ref{doi:10.1093/nar/gkac328}].`)
  .build()
