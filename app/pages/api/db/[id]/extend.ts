import fpprg from '@/app/fpprg'
import type { NextApiRequest, NextApiResponse } from 'next'
import { z } from 'zod'
import { Process } from '@/core/FPPRG'

const QueryType = z.object({
  id: z.string(),
})

const BodyType = z.object({
  type: z.string(),
  inputs: z.record(z.string(), z.string()),
})

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    const { id } = QueryType.parse(req.query)
    const { type, inputs: rawInputs } = BodyType.parse(JSON.parse(req.body))
    const inputs: Record<string, Process> = {}
    for (const arg in rawInputs) {
      const proc = fpprg.getProcess(rawInputs[arg])
      if (proc === undefined) throw new Error('Invalid input')
      inputs[arg] = proc
    }
    const process = fpprg.upsertProcess(new Process(type, inputs))
    const fpl = fpprg.getFPL(id).extend(process)
    res.status(200).json(fpl.toJSON())
  } catch (e) {
    console.error(e)
    res.status(500).end(e.toString())
  }
}
