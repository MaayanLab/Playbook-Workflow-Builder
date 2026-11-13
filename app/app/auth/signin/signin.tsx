'use client'

import React from 'react'
import { useSearchParams } from 'next/navigation'
import { authClient } from '@/app/auth-client'
import { signIn } from '@/app/extensions/auth-client'

export type ProviderDisplay = { name: string, logo?: string }

type Providers = {
  dev: ProviderDisplay | false,
  google: ProviderDisplay | false,
  keycloak: ProviderDisplay | false,
  orcid: ProviderDisplay | false,
  email: boolean,
}

export default function SignIn({ providers }: { providers: Providers }) {
  const searchParams = useSearchParams()
  const callbackURL = searchParams?.get('callbackUrl') ?? '/'
  const [email, setEmail] = React.useState('')
  const [sent, setSent] = React.useState(false)
  const [pending, setPending] = React.useState<string | undefined>()
  const [error, setError] = React.useState<string | undefined>()
  const signInWith = React.useCallback(async (provider: string) => {
    setError(undefined)
    setPending(provider)
    try {
      const result = await signIn(provider, { callbackUrl: callbackURL })
      // (on success the browser is redirected to the provider)
      if (result?.error) setError(result.error.message || result.error.statusText || 'Signing in failed')
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e))
    } finally {
      setPending(undefined)
    }
  }, [callbackURL])
  return (
    <main className="min-h-screen grid place-items-center p-4">
      <div className="card bg-base-100 shadow-xl w-full max-w-md">
        <div className="card-body gap-3">
          <h1 className="card-title">Sign in</h1>
          {(['dev', 'keycloak', 'orcid', 'google'] as const).map(id => {
            const provider = providers[id]
            if (!provider) return null
            return (
              <button key={id} className="btn" disabled={pending !== undefined} onClick={() => signInWith(id)}>
                {provider.logo ? <img src={provider.logo} alt="" className="h-6 w-6 object-contain" /> : null}
                {id === 'dev' ? provider.name : `Sign in with ${provider.name}`}
                {pending === id ? <span className="loading loading-spinner loading-sm" /> : null}
              </button>
            )
          })}
          {providers.email && <form className="flex flex-col gap-2" onSubmit={async event => {
            event.preventDefault()
            setError(undefined)
            const result = await authClient.signIn.magicLink({ email, callbackURL })
            if (result.error) setError(result.error.message || 'Sending the sign-in link failed')
            else setSent(true)
          }}>
            <input className="input input-bordered" type="email" required value={email} onChange={event => setEmail(event.target.value)} placeholder="Email" />
            <button className="btn" type="submit">Email me a sign-in link</button>
            {sent && <p>Check your email for a sign-in link.</p>}
          </form>}
          {error ? <div className="alert alert-error">{error}</div> : null}
        </div>
      </div>
    </main>
  )
}
