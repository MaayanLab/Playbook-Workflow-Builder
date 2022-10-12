import fpprg from '@/app/fpprg'
import type { NextApiRequest, NextApiResponse } from 'next'
import { z } from 'zod'
import { IdOrProcessC } from '@/core/FPPRG'

const QueryType = z.object({
  fpl_id: z.string(),
})

const BodyType = IdOrProcessC

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    if (req.method !== 'POST') throw new Error('Unsupported method')
    const { fpl_id } = QueryType.parse(req.query)
    const process = fpprg.resolveProcess(BodyType.parse(JSON.parse(req.body)))
    const fpl = fpprg.getFPL(fpl_id).extend(process)
    res.status(200).json(fpl.id)
  } catch (e) {
    console.error(e)
    res.status(500).end(e.toString())
  }
}
