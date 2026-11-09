import { getRequestSessionWithId } from '@/app/extensions/auth/server'
import type { SessionWithId } from '@/app/auth'
import type { Context } from './init'

function lazy<T>(fn: () => Promise<T>) {
  let value: Promise<T> | undefined
  return () => value ??= fn()
}

/** Context for server-side callers (e.g. the MCP server) */
export function createServerContext(session: SessionWithId | null = null): Context {
  return { getSession: async () => session }
}

/** Context for App Router Fetch request handlers. */
export function createFetchContext(request: Request): Context {
  return {
    req: { headers: { authorization: request.headers.get('authorization') ?? undefined } },
    getSession: lazy(() => getRequestSessionWithId(request)),
  }
}
