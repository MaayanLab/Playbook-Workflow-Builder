import krg from '@/app/krg'
import fpprg from '@/app/fpprg'
import FPL2BCO from '@/core/fpl2bco'
import { z } from 'zod'
import handler from '@/utils/next-rest'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'

const QueryType = z.object({
  fpl_id: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError()
  const { fpl_id } = QueryType.parse(req.query)
  const fpl = await fpprg.getFPL(fpl_id)
  if (fpl === undefined) throw new NotFoundError()
  const fullFPL = await Promise.all(fpl.resolve().map(fpl => fpl.toJSONWithOutput()))
  const BCO = FPL2BCO(krg, fullFPL)
  res.status(200).json(BCO)
})
