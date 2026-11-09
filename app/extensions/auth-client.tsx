import React from 'react'
import { authClient } from '@/app/auth-client'
import type { SessionWithId } from '@/app/auth'

export function useSession(opts?: { required?: boolean }) {
  const session = authClient.useSession()
  const data: SessionWithId | null = session.data ? {
    user: session.data.user,
    expires: new Date(session.data.session.expiresAt).toISOString(),
  } : null
  React.useEffect(() => {
    if (opts?.required && !session.isPending && !data) {
      window.location.assign(`/auth/signin?callbackUrl=${encodeURIComponent(window.location.href)}`)
    }
  }, [opts?.required, session.isPending, data])
  return { data, status: session.isPending ? 'loading' as const : data ? 'authenticated' as const : 'unauthenticated' as const }
}

export async function signIn(provider?: string, options?: { callbackUrl?: string }) {
  const callbackURL = options?.callbackUrl ?? (typeof window !== 'undefined' ? window.location.href : '/')
  if (!provider) {
    if (typeof window !== 'undefined') window.location.assign(`/auth/signin?callbackUrl=${encodeURIComponent(callbackURL)}`)
    return
  }
  if (provider === 'dev') {
    const result = await authClient.signIn.anonymous()
    if (!result.error && typeof window !== 'undefined') window.location.assign(callbackURL)
    return result
  }
  return await authClient.signIn.social({ provider, callbackURL })
}

export async function signOut(options?: { callbackUrl?: string }) {
  return await authClient.signOut({
    fetchOptions: {
      onSuccess: () => {
        if (typeof window !== 'undefined') window.location.assign(options?.callbackUrl ?? '/')
      },
    },
  })
}

export async function linkSocial(provider: string, callbackURL?: string) {
  return await authClient.linkSocial({
    provider,
    callbackURL: callbackURL ?? (typeof window !== 'undefined' ? window.location.href : '/'),
  })
}

export function SessionProvider({ children }: React.PropsWithChildren<{ session?: unknown }>) {
  return <>{children}</>
}

export type Session = SessionWithId

