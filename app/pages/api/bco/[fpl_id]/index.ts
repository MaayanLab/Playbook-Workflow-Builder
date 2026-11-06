import fpprg from '@/app/fpprg'
import { z } from 'zod'
import handler from '@/utils/next-rest'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import { openApiHandler } from '@/app/openapi'

const QueryType = z.object({
  fpl_id: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'HEAD' && req.method !== 'GET') throw new UnsupportedMethodError(req.method)
  const { fpl_id } = QueryType.parse(req.query)
  if (req.method === 'HEAD') {
    const fpl = await fpprg.getFPL(fpl_id)
    if (fpl === undefined) throw new NotFoundError()
    res.status(200).end()
    return
  }
  req.query.trpc = ['bco', fpl_id]
  await openApiHandler(req, res)
})
