import { z } from 'zod'
import { API } from '@/spec/api'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { UnauthorizedError } from '@/spec/error'
import db from '@/app/db'

export const UserIntegrationsCAVATICA = API('/api/v1/user/integrations/cavatica')
  .query(z.object({}))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    const integrations = await db.objects.user_integrations.findUnique({
      where: {
        id: session.user.id
      }
    })
    if (integrations === null) {
      return {
        id: session.user.id, 
        cavatica_api_key: '',
        cavatica_default_project: '',
      }
    } else {
      return integrations
    }
  })
  .build()

export const UserIntegrationsCAVATICAUpdate = API('/api/v1/user/integrations/cavatica/update')
  .query(z.object({}))
  .body(z.object({
    cavatica_api_key: z.string(),
    cavatica_default_project: z.string(),
  }))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    const integrations = await db.objects.user_integrations.upsert({
      create: {
        id: session.user.id,
        cavatica_api_key: inputs.body.cavatica_api_key,
        cavatica_default_project: inputs.body.cavatica_default_project,
      },
      where: {
        id: session.user.id
      },
      update: {
        cavatica_api_key: inputs.body.cavatica_api_key,
        cavatica_default_project: inputs.body.cavatica_default_project,
      },
    })
    return integrations
  })
  .build()
