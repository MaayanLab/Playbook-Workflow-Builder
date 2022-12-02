import fpprg from '@/app/fpprg'
import { FPL, } from '@/core/FPPRG'
import type { NextApiRequest, NextApiResponse } from 'next'
import { IdOrProcessC } from '@/core/FPPRG'

const BodyType = IdOrProcessC

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    if (req.method !== 'POST') throw new Error('Unsupported method')
    console.log(req.body)
    const process = fpprg.resolveProcess(BodyType.parse(JSON.parse(req.body)))
    const fpl = fpprg.upsertFPL(new FPL(process))
    res.status(200).json(fpl.id)
  } catch (e) {
    console.error(e)
    res.status(500).end((e as Error).toString())
  }
}
