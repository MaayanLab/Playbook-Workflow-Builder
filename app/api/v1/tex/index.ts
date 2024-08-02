import { z } from 'zod'
import { API } from '@/spec/api'
import FPL2TEX from '@/core/fpl2tex'
import fpprg from '@/app/fpprg'
import krg from '@/app/krg'
import { NotFoundError } from '@/spec/error'
import JSZip from 'jszip'
import * as dict from '@/utils/dict'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import db from '@/app/db'

export const TeXForPlaybook = API.get('/api/v1/tex/[fpl_id]')
  .query(z.object({
    fpl_id: z.string(),
    format: z.enum(['zip']).default('zip'),
    metadata: z.string().optional(),
  }))
  .call(async (inputs, req, res) => {
    const fpl_id = inputs.query.fpl_id
    const metadata = z.string().optional().transform(param =>
      param ? z.object({
        title: z.string().optional(),
        description: z.string().optional(),
      }).parse(JSON.parse(param)) : undefined
    ).parse(inputs.query.metadata)
    
    const fpl = await fpprg.getFPL(fpl_id)
    if (fpl === undefined) throw new NotFoundError()
    const session = await getServerSessionWithId(req, res)
    const user = (session && session.user) ? (await db.objects.user.findUnique({ where: { id: session.user.id } })) : undefined
    // @ts-ignore
    const userOrcidAccount = user ? await db.objects.account.findUnique({ where: { userId: user.id, provider: 'orcid' } }) : undefined
    const tex_files = await FPL2TEX({
      krg,
      fpl,
      metadata,
      author: user && user.name !== undefined ? {
        name: user.name,
        ...dict.filter({
          email: user.email,
          affiliation: user.affiliation,
          orcid: userOrcidAccount ? userOrcidAccount.providerAccountId : undefined,
        }),
      } : undefined,
    })
    if (inputs.query.format === 'zip') {
      const zip = new JSZip()
      dict.items(tex_files).forEach(({ key, value }) => {
        zip.file(key, value)
      })
      zip.generateNodeStream().pipe(res)
    }
  })
  .build()
