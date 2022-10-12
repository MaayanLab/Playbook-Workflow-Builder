import fpprg from '@/app/fpprg'
import { FPL, } from '@/core/FPPRG'
import type { NextApiRequest, NextApiResponse } from 'next'
import { z } from 'zod'
import { IdOrProcessC } from '@/core/FPPRG'

const BodyType = z.array(IdOrProcessC)

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    if (req.method !== 'POST') throw new Error('Unsupported method')
    const processArray = BodyType.parse(JSON.parse(req.body)).map(fpprg.resolveProcess)
    const processArrayFPL = FPL.fromProcessArray(processArray)
    if (!processArrayFPL) {
      res.status(404).end()
    } else {
      const fpl = fpprg.upsertFPL(processArrayFPL)
      res.status(200).json(fpl.id)
    }
  } catch (e) {
    console.error(e)
    res.status(500).end((e as Error).toString())
  }
}
