import fpprg from '@/app/fpprg'
import { z } from 'zod'
import handler from '@/utils/next-rest'
import { NotFoundError } from '@/spec/error'

const QueryType = z.object({
  process_id: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'POST') throw new Error('Unsupported method')
  const { process_id } = QueryType.parse(req.query)
  const resolved = await fpprg.getResolved(process_id)
  if (resolved === undefined) throw new NotFoundError()
  await fpprg.deleteResolved(resolved)
  res.status(200).json(null)
})
