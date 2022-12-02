import fpprg from '@/app/fpprg'
import type { NextApiRequest, NextApiResponse } from 'next'
import { z } from 'zod'
import { IdOrProcessC } from '@/core/FPPRG'

const QueryType = z.object({
  fpl_id: z.string(),
  process_id: z.string(),
})

const BodyType = IdOrProcessC

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    if (req.method !== 'POST') throw new Error('Unsupported method')
    const { fpl_id, process_id } = QueryType.parse(req.query)
    const old_process = fpprg.getProcess(process_id)
    if (old_process === undefined) {
        res.status(404).end()
      } else {
      const new_process = fpprg.resolveProcess(BodyType.parse(JSON.parse(req.body)))
      const old_fpl = fpprg.getFPL(fpl_id)
      if (old_fpl === undefined) {
        res.status(404).end()
      } else {
        const new_fpl = old_fpl.rebase(old_process, new_process)
        if (new_fpl === undefined) {
          res.status(404).end()
        } else {
          const fpl = fpprg.upsertFPL(new_fpl)
          res.status(200).json(fpl.id)
        }
      }
    }
  } catch (e) {
    console.error(e)
    res.status(500).end((e as Error).toString())
  }
}
