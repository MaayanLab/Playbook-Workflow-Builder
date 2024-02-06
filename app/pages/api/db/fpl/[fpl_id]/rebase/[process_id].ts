import fpprg from '@/app/fpprg'
import { z } from 'zod'
import { IdOrProcessC } from '@/core/FPPRG'
import handler from '@/utils/next-rest'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'

const QueryType = z.object({
  fpl_id: z.string(),
  process_id: z.string(),
})

const BodyType = IdOrProcessC

export default handler(async (req, res) => {
  if (req.method !== 'POST') throw new UnsupportedMethodError(req.method)
  const { fpl_id, process_id } = QueryType.parse(req.query)
  const old_process = await fpprg.getProcess(process_id)
  if (old_process === undefined) throw new NotFoundError()
  const new_process = await fpprg.resolveProcess(BodyType.parse(req.body))
  const old_fpl = await fpprg.getFPL(fpl_id)
  if (old_fpl === undefined) throw new NotFoundError()
  const { rebased, head } = old_fpl.rebase(old_process, new_process)
  const fpl = await fpprg.upsertFPL(head)
  res.status(200).json({ head: fpl.id, rebased: rebased.id })
})
