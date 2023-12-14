import { z } from 'zod'
import { API } from '@/spec/api'
import db from '@/app/db'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'

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
