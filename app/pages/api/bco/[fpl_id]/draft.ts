import krg from '@/app/krg'
import fpprg from '@/app/fpprg'
import db from '@/app/db'
import FPL2BCO from '@/core/fpl2bco'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import handler from '@/utils/next-rest'
import { NotFoundError, UnsupportedMethodError, UnauthorizedError } from '@/spec/error'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'

const QueryType = z.object({
  fpl_id: z.string(),
  metadata: z.string().optional().transform(param =>
    param ? z.object({
      title: z.string().optional(),
      description: z.string().optional(),
    }).parse(JSON.parse(param)) : undefined
  ),
})

export default handler(async (req, res) => {
  if (req.method !== 'HEAD' && req.method !== 'POST') throw new UnsupportedMethodError(req.method)
  const { fpl_id, metadata } = QueryType.parse(req.query)
  const fpl = await fpprg.getFPL(fpl_id)
  if (fpl === undefined) throw new NotFoundError()
  const session = await getServerSessionWithId(req, res)
  const user = (session && session.user) ? (await db.objects.user.findUnique({ where: { id: session.user.id } })) : undefined
  if (!user) throw new UnauthorizedError()
  if (!user.name) throw new Error('Name required')
  // @ts-ignore
  const userOrcidAccount = await db.objects.account.findUnique({ where: { userId: user.id, provider: 'orcid' } })
  if (!userOrcidAccount) throw new UnauthorizedError('ORCID Required')
  if (Date.now()/1000 > userOrcidAccount.expires_at) throw new UnauthorizedError('ORCID Expired')
  if (req.method === 'HEAD') {
    res.status(200).end()
    return
  }
  const BCO = await FPL2BCO({
    krg,
    fpl,
    metadata,
    author: {
      name: user.name,
      ...dict.filter({
        email: user.email,
        affiliation: user.affiliation,
        orcid: userOrcidAccount.providerAccountId,
      }),
    },
  })
  const bcoReq = await fetch('https://biocomputeobject.org/api/objects/drafts/create/', {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Authorization': `Bearer ${userOrcidAccount.id_token}`,
    },
    body: JSON.stringify({
      POST_api_objects_draft_create: [
        {
          contents: BCO,
          prefix: "BCO",
          schema: "IEEE",
          owner_group: "bco_drafter"
        }
      ]
    }),
  })
  if (bcoReq.status !== 200) {
    res.status(bcoReq.status).end(await bcoReq.text())
    return
  }
  const [bcoRes] = z.array(z.object({ object_id: z.string() })).parse(await bcoReq.json())
  res.status(200).json(bcoRes)
})
