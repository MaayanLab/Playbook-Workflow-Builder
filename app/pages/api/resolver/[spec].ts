import fs from 'fs'
import multiparty from "multiparty"
import type { NextApiRequest, NextApiResponse } from 'next'
import krg from '@/app/krg'
import * as t from 'io-ts'
import decodeOrThrow from '@/utils/decodeOrThrow'

const QueryType = t.type({
  spec: t.string,
})

export const config = {
  api: {
    bodyParser: false,
  },
}

function one<T>(L: T[]): T {
  if (L.length === 0) { throw new Error('Expected value, got none') }
  else if (L.length === 1) return L[0]
  else {
    console.log('Expected one, got multiple, using first')
    return L[0]
  }
}

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    const { spec } = decodeOrThrow(QueryType, req.query)
    const processNode = krg.getResolveNode(spec)
    const form = new multiparty.Form()
    const raw = await new Promise<{ fields: Record<string, string[]>, files: Record<string, multiparty.File[]> }>((resolve, reject) => {
      form.parse(req, function (err, fields, files) {
        if (err) reject({ err })
        else resolve({ fields, files })
      })
    })
    const inputs = {}
    for (const i in processNode.inputs) {
      if (i in raw.fields) {
        inputs[i] = processNode.inputs[i].codec.decode(one(raw.fields[i]))
      } else if (i in raw.files) {
        inputs[i] = processNode.inputs[i].codec.decode(
          await new Promise((resolve, reject) =>
            fs.readFile(one(raw.files[i]).path, (err, data) => {
              if (err) reject(err)
              else resolve(data.toString())
            })
          )
        )
      }
    }
    const ouput = await processNode.resolve({ inputs })
    res.status(200).json(processNode.output.codec.encode(ouput))
  } catch (e) {
    console.error(e)
    res.status(500).end(e.toString())
  }
}