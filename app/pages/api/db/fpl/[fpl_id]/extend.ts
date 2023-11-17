import fpprg from '@/app/fpprg'
import { z } from 'zod'
import { FPL, IdOrProcessC } from '@/core/FPPRG'
import handler from '@/utils/next-rest'
import { NotFoundError, ResponseCodedError, UnsupportedMethodError } from '@/spec/error'

const QueryType = z.object({
  fpl_id: z.string(),
})

const BodyType = IdOrProcessC

export default handler(async (req, res) => {
  if (!fpprg) throw new ResponseCodedError(503, 'Not ready')
  if (req.method !== 'POST') throw new UnsupportedMethodError()
  const { fpl_id } = QueryType.parse(req.query)
  const process = await fpprg.resolveProcess(BodyType.parse(req.body))
  if (fpl_id === 'start') {
    const fpl = await fpprg.upsertFPL(new FPL(process))
    return res.status(200).json(fpl.id)
  }
  const old_fpl = await fpprg.getFPL(fpl_id)
  if (old_fpl === undefined) throw new NotFoundError()
  const fpl = await fpprg.upsertFPL(old_fpl.extend(process))
  res.status(200).json(fpl.id)
})
