import fpprg from '@/app/fpprg'
import { FPL, Process } from '@/core/FPPRG'
import type { NextApiRequest, NextApiResponse } from 'next'
import { z } from 'zod'

const BodyType = z.object({
  type: z.string(),
  inputs: z.record(z.string(), z.string()),
})

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    const { type, inputs: rawInputs } = BodyType.parse(JSON.parse(req.body))
    if (Object.keys(rawInputs).length !== 0) throw new Error('Invalid input')
    const process = fpprg.upsertFPL(new FPL(new Process(type, {})))
    res.status(200).json(process.id)
  } catch (e) {
    console.error(e)
    res.status(500).end(e.toString())
  }
}
