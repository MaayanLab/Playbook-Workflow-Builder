import fpprg from '@/app/fpprg'
import type { NextApiRequest, NextApiResponse } from 'next'
import { z } from 'zod'

const QueryType = z.object({
  process_id: z.string(),
})

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    if (req.method !== 'GET') throw new Error('Unsupported method')
    const { process_id } = QueryType.parse(req.query)
    const process = fpprg.getProcess(process_id)
    if (process === undefined) {
      res.status(404).end()
    } else {
      const output = await process.output()
      res.status(200).json(output ? output.toJSON() : null)
    }
  } catch (e) {
    console.error(e)
    res.status(500).end((e as Error).toString())
  }
}
