import krg from '@/app/krg'
import fpprg from '@/app/fpprg'
import { z } from 'zod'
import handler from '@/utils/next-rest'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import { decode_complete_process_inputs } from '@/core/engine'

const QueryType = z.object({
  process_id: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError(req.method)
  const { process_id } = QueryType.parse(req.query)
  const process = await fpprg.getProcess(process_id)
  if (process === undefined) throw new NotFoundError()
  const inputs = await decode_complete_process_inputs(krg, process)
  res.status(200).json(inputs)
})
