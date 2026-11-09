import SignIn from './signin'
import ClientPage from '@/app/app/_clients/client-page'

export default function SignInPage() {
  return <ClientPage><SignIn providers={{
    dev: process.env.NODE_ENV === 'development',
    google: Boolean(process.env.BETTER_AUTH_GOOGLE ?? process.env.NEXTAUTH_GOOGLE),
    keycloak: Boolean(process.env.BETTER_AUTH_KEYCLOAK ?? process.env.NEXTAUTH_KEYCLOAK),
    orcid: Boolean(process.env.BETTER_AUTH_ORCID ?? process.env.NEXTAUTH_ORCID),
    email: Boolean(process.env.EMAIL_SERVER && process.env.EMAIL_FROM),
  }} /></ClientPage>
}
