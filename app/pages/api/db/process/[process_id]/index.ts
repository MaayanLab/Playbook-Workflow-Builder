import fpprg from '@/app/fpprg'
import { z } from 'zod'
import handler from '@/utils/next-rest'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'

const QueryType = z.object({
  process_id: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError()
  const { process_id } = QueryType.parse(req.query)
  const process = await fpprg.getProcess(process_id)
  if (process === undefined) throw new NotFoundError()
  res.status(200).json(process ? process.toJSON() : null)
})
