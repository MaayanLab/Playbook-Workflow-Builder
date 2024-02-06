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
      const bcoReq = await fetch('https://biocomputeobject.org/users/orcid/user_info/', {
        method: 'POST',
        headers: {
          'Accept': 'application/json',
          'Authorization': `Token ${userOrcidAccount.id_token}`,
          'Referer': 'https://playbook-workflow-builder.cloud',
        },
      })
      return {
        biocompute: bcoReq.status === 200,
        orcid: userOrcidAccount !== undefined,
      }
    }
  })
  .build()

const BCOListResponseC = z.array(z.array(z.object({
  contents: z.object({
    object_id: z.string(),
    provenance_domain: z.object({
      derived_from: z.string(),
    }),
  }),
  state: z.enum(['DRAFT', 'PUBLISHED']),
})))

export const UserIntegrationsBioComputePublishedBCO = API.get('/api/v1/user/integrations/biocompute/[fpl_id]/published')
  .query(z.object({
    fpl_id: z.string(),
  }))
  .call(async (inputs, req, res) => {
    try {
      const session = await getServerSessionWithId(req, res)
      if (!session?.user?.id) throw new UnauthorizedError()
      const user = await db.objects.user.findUnique({ where: { id: session.user.id } })
      if (!user) throw new UnauthorizedError()
      // @ts-ignore
      const userOrcidAccount = user ? await db.objects.account.findUnique({ where: { userId: user.id, provider: 'orcid' } }) : undefined
      if (!userOrcidAccount) throw new UnauthorizedError('ORCID Required')
      // logged in, look for drafts
      const persistent_url = `${process.env.PUBLIC_URL}/report/${inputs.query.fpl_id}`
      const bcoReq = await fetch(`https://biocomputeobject.org/api/objects/?contents=${encodeURIComponent(persistent_url)}`, {
        method: 'GET',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${userOrcidAccount.id_token}`,
        },
      })
      if (!bcoReq.ok) throw new ResponseCodedError(bcoReq.status, await bcoReq.text())
      const results = BCOListResponseC.parse(await bcoReq.json())
        .flatMap(r => r)
        .filter(item => item.contents.provenance_domain.derived_from === persistent_url)
      const published = results.filter(result => result.state === 'PUBLISHED')
      const drafts = results.filter(result => result.state === 'DRAFT')
      return published.length > 0 ? published[0]
        : drafts.length > 0 ? drafts[0]
        : null
    } catch (e) {
      // not logged in, find published BCOs
      const persistent_url = `${process.env.PUBLIC_URL}/report/${inputs.query.fpl_id}`
      const bcoReq = await fetch(`https://biocomputeobject.org/api/objects/?contents=${encodeURIComponent(persistent_url)}`, {
        method: 'GET',
        headers: {
          'Content-Type': 'application/json',
        },
      })
      if (!bcoReq.ok) throw new ResponseCodedError(bcoReq.status, await bcoReq.text())
      const results = BCOListResponseC.parse(await bcoReq.json())
        .flatMap(r => r)
        .filter(item => item.contents.provenance_domain.derived_from === persistent_url)
      return results.length > 0 ? results[0] : null
    }
  })
  .build()
