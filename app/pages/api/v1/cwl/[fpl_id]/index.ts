import { z } from 'zod'
import JSZip from 'jszip'
import YAML from 'yaml'
import { cwl_for_playbook } from '@/core/fpl2cwl'
import fpprg from '@/app/fpprg'
import { getKRG } from '@/app/server-krg'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import handler from '@/utils/next-rest'
import * as dict from '@/utils/dict'
import { openApiHandler } from '@/app/openapi'

const QueryType = z.object({
  fpl_id: z.string(),
  format: z.enum(['zip', 'list', 'json']).default('zip'),
})

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError(req.method)
  const { fpl_id, format } = QueryType.parse(req.query)
  if (format === 'json' || format === 'list') {
    req.query.trpc = ['v1', 'cwl', fpl_id]
    await openApiHandler(req, res)
  } else if (format === 'zip') {
    const fpl = await fpprg.getFPL(fpl_id)
    if (fpl === undefined) throw new NotFoundError()
    const krg = await getKRG()
    const cwl_files = await cwl_for_playbook({ krg, fpl })
    const zip = new JSZip()
    dict.items(cwl_files).forEach(({ key, value }) => {
      zip.file(key, YAML.stringify(value))
    })
    res.setHeader('Content-Type', 'application/zip')
    zip.generateNodeStream().pipe(res)
  }
})
