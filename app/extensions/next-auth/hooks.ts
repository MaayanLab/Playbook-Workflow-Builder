import * as Auth from 'next-auth/react'
import type { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'

export function useSessionWithId(opts: { required: boolean }) {
  const { data, status } = Auth.useSession(opts)
  return { data: data as SessionWithId | null, status }
}
