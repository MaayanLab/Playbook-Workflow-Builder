import { z } from 'zod'
import { API } from '@/spec/api'
import db from '@/app/db'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { ResponseCodedError, UnauthorizedError } from '@/spec/error'

export const UserIntegrationsBioComputeAuth = API.get('/api/v1/user/integrations/biocompute/auth')
  .query(z.object({}))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    const user = (session && session.user) ? (await db.objects.user.findUnique({ where: { id: session.user.id } })) : undefined
    // @ts-ignore
    const userOrcidAccount = user ? await db.objects.account.findUnique({ where: { userId: user.id, provider: 'orcid' } }) : undefined
    if (!userOrcidAccount) {
      return { orcid: false }
    } else {
      // TODO: probably a better endpoint for this
      const bcoReq = await fetch('https://biocomputeobject.org/api/objects/drafts/create/', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${userOrcidAccount.id_token}`,
        },
        body: JSON.stringify({}),
      })
      return {
        biocompute: bcoReq.status !== 403,
        orcid: userOrcidAccount !== undefined,
      }
    }
  })
  .build()

export const UserIntegrationsBioComputePublishedBCO = API.get('/api/v1/user/integrations/biocompute/[fpl_id]/published')
  .query(z.object({
    fpl_id: z.string(),
  }))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session?.user?.id) throw new UnauthorizedError()
    const user = await db.objects.user.findUnique({ where: { id: session.user.id } })
    if (!user) throw new UnauthorizedError()
    // @ts-ignore
    const userOrcidAccount = user ? await db.objects.account.findUnique({ where: { userId: user.id, provider: 'orcid' } }) : undefined
    if (!userOrcidAccount) throw new UnauthorizedError('ORCID Required')
    // TODO: this gives me unauthorized even when ORCID is properly configured
    const bcoReq = await fetch(`https://biocomputeobject.org/api/objects/?contents=${encodeURIComponent(`${process.env.PUBLIC_URL}/report/${inputs.query.fpl_id}`)}`, {
      method: 'GET',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${userOrcidAccount.id_token}`,
      },
    })
    if (!bcoReq.ok) throw new ResponseCodedError(bcoReq.status, await bcoReq.text())
    const results = await bcoReq.json()
    // TODO: verify that the results are as we expect
    // console.log(results)
    return Array.isArray(results) && results.length > 0 ? results[0].object_id : null
  })
  .build()
