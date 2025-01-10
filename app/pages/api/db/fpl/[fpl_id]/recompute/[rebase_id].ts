import fpprg from '@/app/fpprg'
import { z } from 'zod'
import { Process } from '@/core/FPPRG'
import handler from '@/utils/next-rest'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'

const QueryType = z.object({
  fpl_id: z.string(),
  rebase_id: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'POST') throw new UnsupportedMethodError(req.method)
  const { fpl_id, rebase_id } = QueryType.parse(req.query)
  const rebase = await fpprg.getFPL(rebase_id)
  if (rebase === undefined) throw new NotFoundError()
  const old_process = rebase.process
  const old_fpl = await fpprg.getFPL(fpl_id)
  if (old_fpl === undefined) throw new NotFoundError()
  if (old_process.timestamp !== undefined) {
    const new_process = await fpprg.upsertProcess(new Process(
      old_process.type,
      old_process.data,
      old_process.inputs,
    ))
    const { rebased, head } = old_fpl.rebase(old_process, new_process)
    const fpl = await fpprg.upsertFPL(head)
    res.status(200).json({ head: fpl.id, rebased: rebased.id })
  } else {
    const resolved = await fpprg.getResolved(old_process.id)
    if (resolved === undefined) throw new NotFoundError()
    await fpprg.deleteResolved(resolved)
    res.status(200).json({ head: fpl_id, rebased: rebase_id })
  }
})
