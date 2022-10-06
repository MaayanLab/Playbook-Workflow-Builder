import fs from 'fs'
import multiparty from "multiparty"
import type { NextApiRequest, NextApiResponse } from 'next'
import { processNodes } from '@/app/nodes'
import { MetaNode, MetaNodeProcessResolve } from '@/spec/metanode'

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
  const { spec } = req.query
  const processNode = processNodes[spec as string] as MetaNode<MetaNodeProcessResolve>
  const form = new multiparty.Form()
  const raw = await new Promise<{ fields: Record<string, string[]>, files: Record<string, multiparty.File[]> }>((resolve, reject) => {
    form.parse(req, function (err, fields, files) {
      if (err) reject({ err })
      else resolve({ fields, files })
    })
  })
  const inputs = {}
  for (const i in processNode.t.inputs) {
    if (i in raw.fields) {
      inputs[i] = processNode.t.inputs[i].t.codec.decode(one(raw.fields[i]))
    } else if (i in raw.files) {
      inputs[i] = processNode.t.inputs[i].t.codec.decode(
        await new Promise((resolve, reject) =>
          fs.readFile(one(raw.files[i]).path, (err, data) => {
            if (err) reject(err)
            else resolve(data)
          })
        )
      )
    }
  }
  const ouput = await processNode.t.resolve({ inputs })
  res.status(200).json(processNode.t.output.t.codec.encode(ouput))
}