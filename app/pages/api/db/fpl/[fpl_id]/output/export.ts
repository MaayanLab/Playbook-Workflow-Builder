import fpprg from '@/app/fpprg'
import { z } from 'zod'
import handler from '@/utils/next-rest'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import fpl2json from '@/core/fpl2json'

const QueryType = z.object({
  fpl_id: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError(req.method)
  const { fpl_id } = QueryType.parse(req.query)
  const fpl = await fpprg.getFPL(fpl_id)
  if (fpl === undefined) throw new NotFoundError()
  res.status(200).json(await fpl2json({ fpl }))
})
