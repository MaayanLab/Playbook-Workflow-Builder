import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import { ScoredDrugs, ScoredGenes } from '@/components/core/input/scored'

const sigcom_lincs_url = 'https://maayanlab.cloud/sigcom-lincs'

async function sigcom_meta_entities_find_by_symbol(symbols: string[]) {
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
            inq: symbols,
          }
        }
      }
    }),
  })
  return z.array(z.object({
    id: z.string(),
    meta: z.object({
      symbol: z.string(),
    }),
  })).parse(await metaEntityReq.json())
}

async function sigcom_signatures_find_by_id(ids: string[]) {
  const metaSignatureReq = await fetch(`${sigcom_lincs_url}/metadata-api/signatures/find`, {
    method: 'POST',
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'application/json',
    },
    body: JSON.stringify({
      filter: {
        where: {
          id: {
            inq: ids,
          }
        }
      }
    }),
  })
  return z.array(z.object({
    $validator: z.string(),
    id: z.string(),
    meta: z.record(z.string(), z.any()),
  })).parse(await metaSignatureReq.json())
}

async function sigcom_meta_user_input_signature(body: { up_entities: string[], down_entities: string[] }) {
  const metaInputReq = await fetch(`${sigcom_lincs_url}/metadata-api/user_input`, {
    method: 'POST',
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'application/json',
    },
    body: JSON.stringify({
      meta: {
        "$validator": "/dcic/signature-commons-schema/v6/meta/user_input/user_input.json",
        ...body,
      },
      type: "signature",
    }),
  })
  return z.object({ id: z.string() }).parse(await metaInputReq.json())
}

async function sigcom_data_enrich(body: {
  database: string,
  up_entities: string[],
  down_entities: string[],
  limit?: number,
}) {
  const dataReq = await fetch(`${sigcom_lincs_url}/data-api/api/v1/enrich/ranktwosided`, {
    method: 'POST',
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'application/json',
    },
    body: JSON.stringify(body),
  })
  return z.object({
    results: z.array(z.object({
      uuid: z.string(),
      type: z.string(),
      'z-up': z.number(),
      'z-down': z.number(),
    })),
  }).parse(await dataReq.json())
}

export const SigComLINCSSignatureResults = MetaNode(`SigComLINCSSignatureResults`)
  .meta({
    label: 'SigCom LINCS Signature Search Results',
    description: 'LINCS L1000 Signatures Query Results',
  })
  .codec(z.object({ id: z.string(), up_entities: z.array(z.string()), down_entities: z.array(z.string()) }))
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
    const metaRes = await sigcom_meta_entities_find_by_symbol([...dict.keys(up), ...dict.keys(down)])
    const up_entities: string[] = []
    const down_entities: string[] = []
    metaRes.forEach((ent) => {
      if (up[ent.meta.symbol]) up_entities.push(ent.id)
      if (down[ent.meta.symbol]) down_entities.push(ent.id)
    })
    const results = await sigcom_meta_user_input_signature({ up_entities, down_entities })
    return { id: results.id, up_entities, down_entities }
  })
  .story(props => `Identified reversers and mimickers from over 1 million signatures using SigCom LINCS [\\ref{doi:10.1093/nar/gkac328}].`)
  .build()

export const ExtractSigComLINCSSignatureSearchT_l1000_cp = MetaNode(`ExtractSigComLINCSSignatureSearch[LINCS L1000 Chemical Perturbagens]`)
  .meta({
    label: `Extract LINCS L1000 Chemical Perturbagens`,
    description: 'Extract signatures from the results',
  })
  .inputs({ searchResults: SigComLINCSSignatureResults })
  .output(ScoredDrugs)
  .resolve(async (props) => {
    const { results } = await sigcom_data_enrich({
      up_entities: props.inputs.searchResults.up_entities,
      down_entities: props.inputs.searchResults.down_entities,
      limit: 100,
      database: 'l1000_cp',
    })
    const signatureScores = dict.init(results.map(({ uuid: key, ...value }) => ({ key, value })))
    const signatureMeta = await sigcom_signatures_find_by_id(dict.keys(signatureScores))
    const signatures = signatureMeta.map(signature => ({ term: signature.meta.pert_name, zscore: signatureScores[signature.id]['z-up'] - signatureScores[signature.id]['z-down'] }))
    signatures.sort((a, b) => b.zscore - a.zscore)
    return signatures
  })
  .story(props =>
    `Resolved drugs from the LINCS L1000 Chemical Perturbagens library.`
  )
  .build()

export const ExtractSigComLINCSSignatureSearchT_l1000_xpr = MetaNode(`ExtractSigComLINCSSignatureSearch[LINCS L1000 CRISPR KOs]`)
  .meta({
    label: `Extract LINCS L1000 CRISPR KOs`,
    description: 'Extract signatures from the results',
  })
  .inputs({ searchResults: SigComLINCSSignatureResults })
  .output(ScoredGenes)
  .resolve(async (props) => {
    const { results } = await sigcom_data_enrich({
      up_entities: props.inputs.searchResults.up_entities,
      down_entities: props.inputs.searchResults.down_entities,
      limit: 100,
      database: 'l1000_xpr',
    })
    const signatureScores = dict.init(results.map(({ uuid: key, ...value }) => ({ key, value })))
    const signatureMeta = await sigcom_signatures_find_by_id(dict.keys(signatureScores))
    const signatures = signatureMeta.map(signature => ({ term: signature.meta.pert_name, zscore: signatureScores[signature.id]['z-up'] - signatureScores[signature.id]['z-down'] })).filter(({ term }) => !!term)
    signatures.sort((a, b) => b.zscore - a.zscore)
    return signatures
  })
  .story(props =>
    `Resolved genes from the LINCS L1000 CRISPR KOs library.`
  )
  .build()
