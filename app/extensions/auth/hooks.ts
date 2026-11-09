import * as Auth from '@/app/extensions/auth-client'
import type { SessionWithId } from '@/app/auth'

export function useSessionWithId(opts?: { required: boolean }) {
  const { data, status } = Auth.useSession(opts)
  return { data: data as SessionWithId | null, status }
}

