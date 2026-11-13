import SignIn, { type ProviderDisplay } from './signin'
import ClientPage from '@/app/app/_clients/client-page'

/**
 * How a provider configured with a JSON environment variable (e.g. `BETTER_AUTH_KEYCLOAK`) is shown:
 *  its configuration may give a `name` & a logo (`style.logo`), as next-auth's did. Only these are
 *  passed to the browser.
 */
function display(env: string | undefined, defaultName: string): ProviderDisplay | false {
  if (!env) return false
  try {
    const config = JSON.parse(env)
    return {
      name: typeof config.name === 'string' ? config.name : defaultName,
      logo: typeof config.style?.logo === 'string' ? config.style.logo : undefined,
    }
  } catch {
    return { name: defaultName }
  }
}

export default function SignInPage() {
  return <ClientPage><SignIn providers={{
    dev: process.env.NODE_ENV === 'development' ? { name: 'Development user' } : false,
    google: display(process.env.BETTER_AUTH_GOOGLE, 'Google'),
    keycloak: display(process.env.BETTER_AUTH_KEYCLOAK, 'Keycloak'),
    orcid: display(process.env.BETTER_AUTH_ORCID, 'ORCID'),
    email: Boolean(process.env.EMAIL_SERVER && process.env.EMAIL_FROM),
  }} /></ClientPage>
}
