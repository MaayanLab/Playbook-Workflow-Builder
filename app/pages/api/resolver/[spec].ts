import fs from 'fs'
import multiparty from "multiparty"
import type { NextApiRequest, NextApiResponse } from 'next'
import krg from '@/app/krg'
import { z } from 'zod'
import * as array from '@/utils/array'

const QueryType = z.object({
  spec: z.string(),
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
    console.warn('Expected one, got multiple, using first')
    return L[0]
  }
}

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    const { spec } = QueryType.parse(req.query)
    const processNode = krg.getResolveNode(spec)
    const form = new multiparty.Form()
    const raw = await new Promise<{ fields: Record<string, string[]>, files: Record<string, multiparty.File[]> }>((resolve, reject) => {
      form.parse(req, function (err, fields, files) {
        if (err) reject({ err })
        else resolve({ fields, files })
      })
    })
    const inputs: Record<string, unknown> = {}
    for (const i in processNode.inputs) {
      const processNodeInput_i = array.ensureOne(processNode.inputs[i])
      let inputs_i
      if (i in raw.fields) {
        inputs_i = raw.fields[i].map(v => processNodeInput_i.codec.decode(v))
      } else if (i in raw.files) {
        inputs_i = (
          await Promise.all(raw.files[i].map(file =>
            new Promise<string>((resolve, reject) =>
              fs.readFile(file.path, (err, data) => {
                if (err) reject(err)
                else resolve(data.toString())
              })
            )
          ))
        ).map(processNodeInput_i.codec.decode)
      } else {
        throw new Error(`argument ${JSON.stringify(i)} is missing`)
      }
      if (Array.isArray(processNode.inputs[i])) inputs[i] = inputs_i
      else inputs[i] = one(inputs_i)
    }
    const output = await processNode.resolve({ inputs })
    res.status(200).json(processNode.output.codec.encode(output))
  } catch (e) {
    console.error(e)
    res.status(500).json((e as Error).toString())
  }
}