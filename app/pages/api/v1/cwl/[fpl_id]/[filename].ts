import { z } from 'zod'
import YAML from 'yaml'
import { cwl_for_playbook } from '@/core/fpl2cwl'
import fpprg from '@/app/fpprg'
import krg from '@/app/krg'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import handler from '@/utils/next-rest'

const QueryType = z.object({
  fpl_id: z.string(),
  filename: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError(req.method)
  const { fpl_id, filename } = QueryType.parse(req.query)
  const fpl = await fpprg.getFPL(fpl_id)
  if (fpl === undefined) throw new NotFoundError()
  const cwl_files = await cwl_for_playbook({ krg, fpl })
  const cwl_file = cwl_files[filename as keyof typeof cwl_files]
  if (cwl_file === undefined) throw new NotFoundError()
  res.status(200).end(YAML.stringify(cwl_file))
})
