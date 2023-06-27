import { z } from 'zod'
import { API } from '@/spec/api'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { UnauthorizedError, ResponseCodedError, NotFoundError } from '@/spec/error'
import db from '@/app/db'
import emitter from '@/app/emitter'
import { run_wes_worker, abort_wes_worker } from '@/app/extensions/cavatica'

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
    // TODO: spawn with pg-boss
    const proxy_session = await db.objects.proxy_session.create({ data: {} })
    ;(async () => {
      emitter.on(`join:${proxy_session.id}`, async () => {
        await db.objects.proxy_session.update({
          where: { id: proxy_session.id },
          data: { state: 'CONNECTED' },
        })
      })
      try {
        for await (const status of run_wes_worker({
          auth_token: integrations.cavatica_api_key,
          project: integrations.cavatica_default_project,
          socket: process.env.PUBLIC_URL as string,
          session_id: proxy_session.id,
        })) {
          await db.objects.proxy_session.update({
            where: { id: proxy_session.id },
            data: status,
          })
        }
      } catch (e) {
        await db.objects.proxy_session.update({
          where: { id: proxy_session.id },
          data: { state: 'ERROR' },
        })
      }
      await db.objects.proxy_session.delete({ where: { id: proxy_session.id } })
    })()
    return proxy_session.id
  })
  .build()

export const UserIntegrationsCAVATICAStatus = API('/api/v1/user/integrations/cavatica/[session_id]/status')
  .query(z.object({
    session_id: z.string(),
  }))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    const proxy_session = await db.objects.proxy_session.findUnique({ where: { id: inputs.query.session_id } })
    if (proxy_session === null) throw new NotFoundError()
    return proxy_session
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
    const proxy_session = await db.objects.proxy_session.findUnique({ where: { id: inputs.query.session_id } })
    if (proxy_session === null) throw new NotFoundError()
    if (!proxy_session.run_id) throw new ResponseCodedError(402, 'Not launched yet!')
    if (proxy_session.state === 'RUNNING') {
      emitter.emit(`close:${inputs.query.session_id}`)
    }
    const integrations = await db.objects.user_integrations.findUnique({
      where: {
        id: session.user.id
      }
    })
    if (!integrations?.cavatica_api_key) throw new ResponseCodedError(402, 'CAVATICA Integration not configured')
    await abort_wes_worker({
      run_id: proxy_session.run_id,
      auth_token: integrations.cavatica_api_key,
    })
    return null
  })
  .build()
