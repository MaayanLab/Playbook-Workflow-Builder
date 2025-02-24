import { z } from 'zod'
import { API } from '@/spec/api'
import { fpl2ro_crate, fpl2ro_crate_metadata } from '@/core/fpl2ro-crate'
import fpprg from '@/app/fpprg'
import krg from '@/app/krg'
import { NotFoundError } from '@/spec/error'
import JSZip from 'jszip'
import * as dict from '@/utils/dict'

export const ROCrateForPlaybook = API.get('/api/v1/ro-crate/[fpl_id]')
  .query(z.object({
    fpl_id: z.string(),
    format: z.enum(['zip', 'json']).default('zip'),
  }))
  .call(async (inputs, req, res) => {
    const fpl = await fpprg.getFPL(inputs.query.fpl_id)
    if (fpl === undefined) throw new NotFoundError()
    if (inputs.query.format === 'json') {
      return await fpl2ro_crate_metadata({ krg, fpl })
    } else if (inputs.query.format === 'zip') {
      const zip = new JSZip()
      const ro_crate_files = await fpl2ro_crate({ krg, fpl })
      dict.items(ro_crate_files).forEach(({ key, value }) => {
        zip.file(key, value)
      })
      zip.generateNodeStream().pipe(res)
    }
  })
  .build()
