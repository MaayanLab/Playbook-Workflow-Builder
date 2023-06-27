import { z } from 'zod'
import { API } from '@/spec/api'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { UnauthorizedError, ResponseCodedError } from '@/spec/error'
import { v4 as uuidv4 } from 'uuid'
import db from '@/app/db'
import emitter from '@/app/emitter'
import { run_wes_worker } from '@/app/extensions/cavatica'

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

export const UserIntegrationsCAVATICALaunch = API('/api/v1/user/integrations/cavatica/launch')
  .query(z.object({}))
  .body(z.object({}))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    const integrations = await db.objects.user_integrations.findUnique({
      where: { id: session.user.id }
    })
    if (!integrations?.cavatica_api_key || !integrations?.cavatica_default_project) throw new ResponseCodedError(402, 'CAVATICA Integration not configured')
    const session_id = uuidv4()
    // TODO: spawn with pg-boss
    ;(async () => {
      for await (const status of run_wes_worker({
        auth_token: integrations.cavatica_api_key,
        project: integrations.cavatica_default_project,
        socket: process.env.PUBLIC_URL as string,
        session_id,
      })) {
        if (status.state === null) {
          console.debug(`Started task for user ${session.user.id} with run_id=${status.run_id}`)
        } else {
          console.debug(`Task for user ${session.user.id} entered state: ${status.state}`)
        }
      }
    })()
    return session_id
  })
  .build()

export const UserIntegrationsCAVATICADisconnect = API('/api/v1/user/integrations/cavatica/[session_id]/disconnect')
  .query(z.object({
    session_id: z.string(),
  }))
  .body(z.object({}))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    emitter.emit(`close:${inputs.query.session_id}`)
    return null
  })
  .build()
