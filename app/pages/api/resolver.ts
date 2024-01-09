import type { NextApiRequest, NextApiResponse } from 'next'
import krg from '@/app/krg'
import { z } from 'zod'

const BodyType = z.object({
  spec: z.string(),
  data: z.unknown().optional(),
  inputs: z.record(z.string(), z.unknown()),
})

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    const { spec, data, inputs } = BodyType.parse(req.body)
    const processNode = krg.getProcessNode(spec)
    const output = await processNode.resolve({ data, inputs, notify: (update) => {} })
    res.status(200).json(processNode.output.codec.encode(output))
  } catch (e) {
    console.error(e)
    res.status(500).json((e as Error).toString())
  }
}
