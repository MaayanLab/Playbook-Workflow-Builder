import { z } from 'zod'
import { TRPCError } from '@trpc/server'
import db from '@/app/db'
import { uuid } from '@/db'
import { closeSession, isWorkerConnected, sessionEvents } from '@/app/extensions/cavatica/tunnel'
import * as cavatica from '@/app/extensions/cavatica'
import * as cavaticaLocal from '@/app/extensions/cavatica/testing'
import { getUser, getUserOrcidAccount } from '@/app/extensions/next-auth/user'
import krg from '@/app/krg'
import fpprg from '@/app/fpprg'
import FPL2BCO from '@/core/fpl2bco'
import * as dict from '@/utils/dict'
import { authedProcedure, publicProcedure, router } from '../init'

async function getIntegrations(user_id: string) {
  return await db.selectFrom('user_integrations').selectAll().where('id', '=', user_id).executeTakeFirst()
}

/** A CAVATICA session launched by `user` */
async function getOwnSession(session_id: string, user: string) {
  const proxy_session = await db.selectFrom('proxy_session').selectAll()
    .where('id', '=', session_id)
    .where('user', '=', user)
    .executeTakeFirst()
  if (!proxy_session) throw new TRPCError({ code: 'NOT_FOUND' })
  return proxy_session
}

/** CAVATICA_LOCAL_WORKER=1 runs session workers on this machine instead of in CAVATICA (for development) */
const { run_wes_worker, abort_wes_worker } = process.env.CAVATICA_LOCAL_WORKER ? cavaticaLocal : cavatica

const notConfigured = () => new TRPCError({ code: 'PAYMENT_REQUIRED', message: 'CAVATICA Integration not configured' })

const cavaticaRouter = router({
  get: authedProcedure
    .query(async ({ ctx }) => {
      const integrations = await getIntegrations(ctx.session.user.id)
      return {
        id: ctx.session.user.id,
        cavatica_api_key: integrations?.cavatica_api_key ?? '',
        cavatica_default_project: integrations?.cavatica_default_project ?? '',
      }
    }),

  update: authedProcedure
    .input(z.object({
      cavatica_api_key: z.string(),
      cavatica_default_project: z.string(),
    }))
    .mutation(async ({ ctx, input }) => {
      return await db.insertInto('user_integrations')
        .values({ id: ctx.session.user.id, ...input })
        .onConflict(oc => oc.column('id').doUpdateSet(input))
        .returningAll()
        .executeTakeFirstOrThrow()
    }),

  launch: authedProcedure
    .mutation(async ({ ctx }) => {
      const integrations = await getIntegrations(ctx.session.user.id)
      if (!integrations?.cavatica_api_key || !integrations?.cavatica_default_project) throw notConfigured()
      const { cavatica_api_key, cavatica_default_project } = integrations
      // TODO: spawn with pg-boss
      const proxy_session = await db.insertInto('proxy_session')
        .values({ id: uuid(), user: ctx.session.user.id, state: '' })
        .returning('id')
        .executeTakeFirstOrThrow()
      const setState = (data: { state?: string | null, run_id?: string | null }) =>
        db.updateTable('proxy_session').set(data).where('id', '=', proxy_session.id).execute()
      ;(async () => {
        const onJoin = () => { setState({ state: 'CONNECTED' }).catch(error => console.error(error)) }
        sessionEvents.on(`join:${proxy_session.id}`, onJoin)
        try {
          for await (const status of run_wes_worker({
            url: process.env.PUBLIC_URL as string,
            auth_token: cavatica_api_key,
            project: cavatica_default_project,
            session_id: proxy_session.id,
          })) {
            await setState(status)
          }
        } catch (e) {
          await setState({ state: 'ERROR' })
        }
        sessionEvents.off(`join:${proxy_session.id}`, onJoin)
        await db.deleteFrom('proxy_session').where('id', '=', proxy_session.id).execute()
      })().catch(error => console.error(error))
      return proxy_session.id
    }),

  status: authedProcedure
    .input(z.object({ session_id: z.string() }))
    .query(async ({ ctx, input }) => {
      return await getOwnSession(input.session_id, ctx.session.user.id)
    }),

  disconnect: authedProcedure
    .input(z.object({ session_id: z.string() }))
    .mutation(async ({ ctx, input }) => {
      const proxy_session = await getOwnSession(input.session_id, ctx.session.user.id)
      if (!proxy_session.run_id) throw new TRPCError({ code: 'PAYMENT_REQUIRED', message: 'Not launched yet!' })
      if (isWorkerConnected(input.session_id)) closeSession(input.session_id)
      const integrations = await getIntegrations(ctx.session.user.id)
      if (!integrations?.cavatica_api_key) throw notConfigured()
      await abort_wes_worker({
        run_id: proxy_session.run_id,
        auth_token: integrations.cavatica_api_key,
      })
      return null
    }),
})

const BCOListResponseC = z.array(z.array(z.object({
  contents: z.object({
    object_id: z.string(),
    provenance_domain: z.object({
      derived_from: z.string(),
    }),
  }),
  state: z.enum(['DRAFT', 'PUBLISHED']),
})))

async function fetchBCOs(fpl_id: string, id_token?: string | null) {
  const persistent_url = `${process.env.PUBLIC_URL}/report/${fpl_id}`
  const bcoReq = await fetch(`https://biocomputeobject.org/api/objects/?contents=${encodeURIComponent(persistent_url)}`, {
    method: 'GET',
    headers: {
      'Content-Type': 'application/json',
      ...(id_token ? { 'Authorization': `Bearer ${id_token}` } : {}),
    },
  })
  if (!bcoReq.ok) throw new TRPCError({ code: 'BAD_GATEWAY', message: await bcoReq.text() })
  return BCOListResponseC.parse(await bcoReq.json())
    .flatMap(r => r)
    .filter(item => item.contents.provenance_domain.derived_from === persistent_url)
}

const biocomputeRouter = router({
  draft: authedProcedure
    .input(z.object({
      fpl_id: z.string(),
      metadata: z.object({ title: z.string().optional(), description: z.string().optional() }).optional(),
    }))
    .mutation(async ({ ctx, input }) => {
      const [fpl, user, account] = await Promise.all([
        fpprg.getFPL(input.fpl_id),
        getUser(ctx.session.user.id),
        getUserOrcidAccount(ctx.session.user.id),
      ])
      if (!fpl) throw new TRPCError({ code: 'NOT_FOUND' })
      if (!user?.name) throw new TRPCError({ code: 'BAD_REQUEST', message: 'Name required' })
      if (!account) throw new TRPCError({ code: 'UNAUTHORIZED', message: 'ORCID Required' })
      if (Date.now() / 1000 > (account.expires_at ?? 0)) throw new TRPCError({ code: 'UNAUTHORIZED', message: 'ORCID Expired' })
      const contents = await FPL2BCO({
        krg,
        fpl,
        metadata: input.metadata,
        author: {
          name: user.name,
          ...dict.filter({ email: user.email, affiliation: user.affiliation ?? undefined, orcid: account.providerAccountId }),
        },
      })
      const response = await fetch('https://biocomputeobject.org/users/bcodb/draft_bco/add', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json', Authorization: `Bearer ${account.id_token}` },
        body: JSON.stringify({ contents, origin: 'https://playbook-workflow-builder.cloud' }),
      })
      if (response.status === 403) throw new TRPCError({ code: 'UNAUTHORIZED', message: 'BCO Unauthorization' })
      if (!response.ok) throw new TRPCError({ code: 'BAD_GATEWAY', message: await response.text() })
      return z.object({ object_id: z.string() }).parse(await response.json())
    }),

  auth: publicProcedure
    .query(async ({ ctx }) => {
      const session = await ctx.getSession()
      const user = session?.user ? await getUser(session.user.id) : undefined
      const userOrcidAccount = user ? await getUserOrcidAccount(user.id) : undefined
      if (!userOrcidAccount) return { orcid: false, biocompute: false }
      const bcoReq = await fetch('https://biocomputeobject.org/users/orcid/user_info/', {
        method: 'POST',
        headers: {
          'Accept': 'application/json',
          'Authorization': `Token ${userOrcidAccount.id_token}`,
          'Referer': 'https://playbook-workflow-builder.cloud',
        },
      })
      return { orcid: true, biocompute: bcoReq.status === 200 }
    }),

  /** The published BCO for this playbook, or the user's draft if there's no published BCO */
  published: publicProcedure
    .input(z.object({ fpl_id: z.string() }))
    .query(async ({ ctx, input }) => {
      const session = await ctx.getSession()
      const userOrcidAccount = session?.user ? await getUserOrcidAccount(session.user.id) : undefined
      if (userOrcidAccount) {
        try {
          // logged in, look for drafts
          const results = await fetchBCOs(input.fpl_id, userOrcidAccount.id_token)
          return results.find(result => result.state === 'PUBLISHED')
            ?? results.find(result => result.state === 'DRAFT')
            ?? null
        } catch (e) {
          console.warn(e)
        }
      }
      // not logged in, find published BCOs
      const results = await fetchBCOs(input.fpl_id)
      return results[0] ?? null
    }),
})

export const integrationsRouter = router({
  cavatica: cavaticaRouter,
  biocompute: biocomputeRouter,
})
