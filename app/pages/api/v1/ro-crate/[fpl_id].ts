import { z } from 'zod'
import JSZip from 'jszip'
import { fpl2ro_crate, fpl2ro_crate_metadata } from '@/core/fpl2ro-crate'
import fpprg from '@/app/fpprg'
import krg from '@/app/krg'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import handler from '@/utils/next-rest'
import * as dict from '@/utils/dict'

const QueryType = z.object({
  fpl_id: z.string(),
  format: z.enum(['zip', 'json']).default('zip'),
})

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError(req.method)
  const { fpl_id, format } = QueryType.parse(req.query)
  const fpl = await fpprg.getFPL(fpl_id)
  if (fpl === undefined) throw new NotFoundError()
  if (format === 'json') {
    res.status(200).json(await fpl2ro_crate_metadata({ krg, fpl }))
  } else if (format === 'zip') {
    const zip = new JSZip()
    const ro_crate_files = await fpl2ro_crate({ krg, fpl })
    dict.items(ro_crate_files).forEach(({ key, value }) => {
      zip.file(key, value)
    })
    res.setHeader('Content-Type', 'application/zip')
    zip.generateNodeStream().pipe(res)
  }
})
