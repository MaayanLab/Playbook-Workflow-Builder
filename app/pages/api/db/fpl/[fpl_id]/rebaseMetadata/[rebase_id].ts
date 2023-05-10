import fpprg from '@/app/fpprg'
import { z } from 'zod'
import { IdOrFPLMetadataC } from '@/core/FPPRG'
import handler from '@/utils/next-rest'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'

const QueryType = z.object({
  fpl_id: z.string(),
  rebase_id: z.string(),
})

const BodyType = IdOrFPLMetadataC

export default handler(async (req, res) => {
  if (req.method !== 'POST') throw new UnsupportedMethodError()
  const { fpl_id, rebase_id } = QueryType.parse(req.query)
  const rebase_fpl = await fpprg.getFPL(rebase_id)
  if (rebase_fpl === undefined) throw new NotFoundError()
  const new_metadata = await fpprg.resolveFPLMetadata(BodyType.parse(JSON.parse(req.body)))
  const old_fpl = await fpprg.getFPL(fpl_id)
  if (old_fpl === undefined) throw new NotFoundError()
  const { rebased, head } = old_fpl.rebaseMetadata(rebase_fpl, new_metadata)
  const fpl = await fpprg.upsertFPL(head)
  res.status(200).json({ head: fpl.id, rebased: rebased.id })
})
