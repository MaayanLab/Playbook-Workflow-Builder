import type { CreateNextContextOptions } from '@trpc/server/adapters/next'
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

/** Context for server-side callers (e.g. the MCP server) */
export function createServerContext(session: SessionWithId | null = null): Context {
  return { getSession: async () => session }
}
