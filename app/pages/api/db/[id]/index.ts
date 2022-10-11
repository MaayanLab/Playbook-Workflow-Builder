import fpprg from '@/app/fpprg'
import type { NextApiRequest, NextApiResponse } from 'next'
import { z } from 'zod'

const QueryType = z.object({
  id: z.string(),
})

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    const { id } = QueryType.parse(req.query)
    const fpl = fpprg.getFPL(id)
    if (fpl === undefined) {
      res.status(404).end()
    } else {
      res.status(200).json(fpl.resolve().map(fpl => fpl.toJSON()))
    }
  } catch (e) {
    console.error(e)
    res.status(500).end(e.toString())
  }
}
