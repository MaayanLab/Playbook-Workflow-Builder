import { initTRPC, TRPCError } from '@trpc/server'
import type { TRPC_ERROR_CODE_KEY } from '@trpc/server/rpc'
import type { IncomingMessage } from 'http'
import { ResponseCodedError } from '@/spec/error'
import { SYSTEM_USER_ID } from '@/db/constants'
import type { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'
import type { OpenApiMeta } from 'trpc-to-openapi'

export type Context = {
  req?: IncomingMessage
  /** resolved lazily, only procedures which need the session pay for it */
  getSession: () => Promise<SessionWithId | null>
}

const t = initTRPC.context<Context>().meta<OpenApiMeta>().create()

const codeForStatus: Record<number, TRPC_ERROR_CODE_KEY> = {
  400: 'BAD_REQUEST',
  401: 'UNAUTHORIZED',
  402: 'PAYMENT_REQUIRED',
  403: 'FORBIDDEN',
  404: 'NOT_FOUND',
  405: 'METHOD_NOT_SUPPORTED',
  408: 'TIMEOUT',
  409: 'CONFLICT',
  412: 'PRECONDITION_FAILED',
  413: 'PAYLOAD_TOO_LARGE',
  415: 'UNSUPPORTED_MEDIA_TYPE',
  422: 'UNPROCESSABLE_CONTENT',
  429: 'TOO_MANY_REQUESTS',
  501: 'NOT_IMPLEMENTED',
  502: 'BAD_GATEWAY',
  503: 'SERVICE_UNAVAILABLE',
  504: 'GATEWAY_TIMEOUT',
}

/**
 * Translate our ResponseCodedErrors (NotFoundError, UnauthorizedError, ...)
 *  into TRPCErrors with the corresponding code (and so http status).
 */
const mapErrors = t.middleware(async ({ next }) => {
  const result = await next()
  if (!result.ok && ResponseCodedError.isinstance(result.error.cause)) {
    const cause = result.error.cause as unknown as ResponseCodedError
    throw new TRPCError({
      code: codeForStatus[cause.error_code] ?? 'INTERNAL_SERVER_ERROR',
      message: cause.message,
      cause,
    })
  }
  return result
})

export const router = t.router
export const createCallerFactory = t.createCallerFactory

export const publicProcedure = t.procedure.use(mapErrors)

export const authedProcedure = publicProcedure.use(async ({ ctx, next }) => {
  const session = await ctx.getSession()
  if (!session?.user) throw new TRPCError({ code: 'UNAUTHORIZED' })
  return next({ ctx: { session: session as SessionWithId & { user: NonNullable<SessionWithId['user']> } } })
})

export const adminProcedure = authedProcedure.use(async ({ ctx, next }) => {
  if (ctx.session.user.id !== SYSTEM_USER_ID) throw new TRPCError({ code: 'UNAUTHORIZED' })
  return next()
})
