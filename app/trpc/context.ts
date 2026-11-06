import type { CreateNextContextOptions } from '@trpc/server/adapters/next'
import type { CreateWSSContextFnOptions } from '@trpc/server/adapters/ws'
import type { IncomingMessage } from 'http'
import { getToken } from 'next-auth/jwt'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import type { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'
import type { Context } from './init'

function lazy<T>(fn: () => Promise<T>) {
  let value: Promise<T> | undefined
  return () => value ??= fn()
}

/** Context for http requests (through the next.js api route) */
export function createNextContext({ req, res }: CreateNextContextOptions): Context {
  return { req, getSession: lazy(() => getServerSessionWithId(req, res)) }
}

/**
 * Resolve the next-auth session from the cookies on a websocket upgrade request.
 *  (we use the jwt session strategy so this doesn't require a response object)
 */
export async function getSessionFromRequest(req: IncomingMessage): Promise<SessionWithId | null> {
  if (process.env.NEXTAUTH_SECRET && req.headers.authorization === `Token ${process.env.NEXTAUTH_SECRET}`) {
    return {
      user: { id: '00000000-0000-0000-0000-000000000000' },
      expires: (new Date(Date.now() + 1000*60*60*24)).toISOString(),
    } as SessionWithId
  }
  const token = await getToken({ req: req as Parameters<typeof getToken>[0]['req'], secret: process.env.NEXTAUTH_SECRET })
  if (!token?.sub || !token.exp) return null
  return {
    user: { id: token.sub, name: token.name, email: token.email, image: token.picture },
    expires: new Date((token.exp as number) * 1000).toISOString(),
  }
}

/** Context for websocket connections */
export function createWSContext({ req }: CreateWSSContextFnOptions): Context {
  return { req, getSession: lazy(() => getSessionFromRequest(req)) }
}

/** Context for server-side callers (e.g. the MCP server) */
export function createServerContext(session: SessionWithId | null = null): Context {
  return { getSession: async () => session }
}
