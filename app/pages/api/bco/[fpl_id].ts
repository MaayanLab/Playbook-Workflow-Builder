import krg from '@/app/krg'
import fpprg from '@/app/fpprg'
import db from '@/app/db'
import FPL2BCO from '@/core/fpl2bco'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import handler from '@/utils/next-rest'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'

const QueryType = z.object({
  fpl_id: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError()
  const { fpl_id } = QueryType.parse(req.query)
  const fpl = await fpprg.getFPL(fpl_id)
  if (fpl === undefined) throw new NotFoundError()
  const fullFPL = await Promise.all(fpl.resolve().map((fpl) => fpl.toJSONWithOutput()))
  const session = await getServerSessionWithId(req, res)
  const user = (session && session.user) ? (await db.objects.user.findUnique({ where: { id: session.user.id } })) : undefined
  const BCO = FPL2BCO({
    krg,
    fpl: fullFPL,
    author: user && user.name !== undefined ? {
      name: user.name,
      ...dict.filter({
        email: user.email,
        affiliation: user.affiliation,
      }),
    } : undefined,
  })
  res.status(200).json(BCO)
})
