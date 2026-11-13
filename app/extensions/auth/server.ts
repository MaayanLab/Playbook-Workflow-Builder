import { SYSTEM_USER_ID } from '@/db/constants'
import { auth, type SessionWithId } from '@/app/auth'

export async function getRequestSessionWithId(request: Request): Promise<SessionWithId | null> {
  const workerToken = process.env.WORKER_AUTH_TOKEN ?? process.env.BETTER_AUTH_SECRET
  if (workerToken && request.headers.get('authorization') === `Token ${workerToken}`) {
    return {
      user: { id: SYSTEM_USER_ID },
      expires: new Date(Date.now() + 1000*60*60*24).toISOString(),
    } as SessionWithId
  }
  // (an invalid or expired API key is an error, rather than no session)
  const result = await auth.api.getSession({ headers: request.headers }).catch(() => null)
  if (!result) return null
  return {
    user: result.user,
    expires: new Date(result.session.expiresAt).toISOString(),
  } as SessionWithId
}
