import fpprg from '@/app/fpprg'
import type { NextApiRequest, NextApiResponse } from 'next'
import { z } from 'zod'
import { FPL, IdOrProcessC } from '@/core/FPPRG'

const QueryType = z.object({
  fpl_id: z.string(),
})

const BodyType = IdOrProcessC

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    if (!fpprg) throw new Error('Not ready')
    if (req.method !== 'POST') throw new Error('Unsupported method')
    const { fpl_id } = QueryType.parse(req.query)
    const process = fpprg.resolveProcess(BodyType.parse(JSON.parse(req.body)))
    if (fpl_id === 'start') {
      const fpl = fpprg.upsertFPL(new FPL(process))
      res.status(200).json(fpl.id)
    } else {
      const old_fpl = fpprg.getFPL(fpl_id)
      if (old_fpl === undefined) {
        res.status(404).end()
      } else {
        const fpl = fpprg.upsertFPL(old_fpl.extend(process))
        res.status(200).json(fpl.id)
      }
    }
  } catch (e) {
    console.error(e)
    res.status(500).end((e as Error).toString())
  }
}
