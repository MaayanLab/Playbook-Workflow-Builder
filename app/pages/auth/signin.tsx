import React from 'react'
import { useRouter } from 'next/router'
import { authClient } from '@/app/auth-client'
import { signIn } from '@/app/extensions/auth-client'

export async function getServerSideProps() {
  return { props: {
    providers: {
      dev: process.env.NODE_ENV === 'development',
      google: Boolean(process.env.BETTER_AUTH_GOOGLE ?? process.env.NEXTAUTH_GOOGLE),
      keycloak: Boolean(process.env.BETTER_AUTH_KEYCLOAK ?? process.env.NEXTAUTH_KEYCLOAK),
      orcid: Boolean(process.env.BETTER_AUTH_ORCID ?? process.env.NEXTAUTH_ORCID),
      email: Boolean(process.env.EMAIL_SERVER && process.env.EMAIL_FROM),
    },
  } }
}

export default function SignIn({ providers }: { providers: Record<string, boolean> }) {
  const router = useRouter()
  const callbackURL = typeof router.query.callbackUrl === 'string' ? router.query.callbackUrl : '/'
  const [email, setEmail] = React.useState('')
  const [sent, setSent] = React.useState(false)
  return (
    <main className="min-h-screen grid place-items-center p-4">
      <div className="card bg-base-100 shadow-xl w-full max-w-md">
        <div className="card-body gap-3">
          <h1 className="card-title">Sign in</h1>
          {providers.dev && <button className="btn" onClick={() => signIn('dev', { callbackUrl: callbackURL })}>Development user</button>}
          {providers.google && <button className="btn" onClick={() => signIn('google', { callbackUrl: callbackURL })}>Continue with Google</button>}
          {providers.keycloak && <button className="btn" onClick={() => signIn('keycloak', { callbackUrl: callbackURL })}>Continue with Keycloak</button>}
          {providers.orcid && <button className="btn" onClick={() => signIn('orcid', { callbackUrl: callbackURL })}>Continue with ORCID</button>}
          {providers.email && <form className="flex flex-col gap-2" onSubmit={async event => {
            event.preventDefault()
            const result = await authClient.signIn.magicLink({ email, callbackURL })
            if (!result.error) setSent(true)
          }}>
            <input className="input input-bordered" type="email" required value={email} onChange={event => setEmail(event.target.value)} placeholder="Email" />
            <button className="btn" type="submit">Email me a sign-in link</button>
            {sent && <p>Check your email for a sign-in link.</p>}
          </form>}
        </div>
      </div>
    </main>
  )
}
