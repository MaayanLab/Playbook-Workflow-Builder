import fpprg from '@/app/fpprg'
import { DataC, IdOrPlaybookMetadataC } from '@/core/FPPRG'
import { z } from 'zod'
import { IdOrProcessC } from '@/core/FPPRG'
import { UnsupportedMethodError } from '@/spec/error'
import handler from '@/utils/next-rest'

const BodyType = z.object({
  data: z.record(z.string(), DataC).optional(),
  workflow: z.array(IdOrProcessC),
  metadata: IdOrPlaybookMetadataC.optional(),
})

export default handler(async (req, res) => {
  if (req.method !== 'POST') throw new UnsupportedMethodError(req.method)
  const { data, workflow, metadata } = BodyType.parse(req.body)
  const fpl = await fpprg.resolveFPL({ data, workflow, metadata })
  res.status(200).json(fpl.id)
})
