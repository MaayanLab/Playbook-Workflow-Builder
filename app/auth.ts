import { betterAuth } from 'better-auth'
import { anonymous, genericOAuth, keycloak, magicLink } from 'better-auth/plugins'
import { apiKey } from '@better-auth/api-key'
import nodemailer from 'nodemailer'
import db, { dialect } from '@/app/db'

type ProviderOptions = { clientId: string, clientSecret: string, issuer?: string }
const parseProvider = (value?: string): ProviderOptions | undefined => value ? JSON.parse(value) : undefined

const google = parseProvider(process.env.BETTER_AUTH_GOOGLE)
const keycloakOptions = parseProvider(process.env.BETTER_AUTH_KEYCLOAK)
const orcid = parseProvider(process.env.BETTER_AUTH_ORCID)
const oauthProviders = [
  // (the configuration may also have display options, e.g. `name` & `style`, see app/app/auth/signin/page.tsx)
  keycloakOptions?.issuer ? keycloak({ clientId: keycloakOptions.clientId, clientSecret: keycloakOptions.clientSecret, issuer: keycloakOptions.issuer }) : undefined,
  orcid ? {
    providerId: 'orcid' as const,
    name: 'ORCID',
    clientId: orcid.clientId,
    clientSecret: orcid.clientSecret,
    discoveryUrl: 'https://orcid.org/.well-known/openid-configuration',
    requireIdTokenVerification: true,
    scopes: ['openid'],
    mapProfileToUser: (profile: Record<string, unknown>) => ({
      name: [profile.given_name, profile.family_name].filter(Boolean).join(' ') || String(profile.sub),
      email: typeof profile.email === 'string' ? profile.email : `${profile.sub}@orcid.invalid`,
    }),
  } : undefined,
].filter((provider): provider is NonNullable<typeof provider> => provider !== undefined)

const emailTransport = process.env.EMAIL_SERVER ? nodemailer.createTransport(process.env.EMAIL_SERVER) : undefined

export const auth = betterAuth({
  // Better Auth's built-in adapter is Kysely-based. Supplying the explicit
  // dialect avoids relying on runtime detection of an already-created Kysely.
  database: { db, type: dialect },
  secret: process.env.BETTER_AUTH_SECRET,
  baseURL: process.env.BETTER_AUTH_URL ?? process.env.PUBLIC_URL,
  user: {
    fields: {
      emailVerified: 'authEmailVerified',
      createdAt: 'created',
      updatedAt: 'authUpdated',
    },
  },
  session: { modelName: 'auth_session' },
  account: { modelName: 'auth_account' },
  verification: { modelName: 'auth_verification' },
  socialProviders: google ? { google } : {},
  plugins: [
    ...(oauthProviders.length > 0 ? [genericOAuth({ config: oauthProviders })] : []),
    ...(emailTransport && process.env.EMAIL_FROM ? [magicLink({
      sendMagicLink: async ({ email, url }) => {
        await emailTransport.sendMail({
          from: process.env.EMAIL_FROM,
          to: email,
          subject: 'Sign in to Playbook Workflow Builder',
          text: `Sign in: ${url}`,
          html: `<p><a href="${url}">Sign in to Playbook Workflow Builder</a></p>`,
        })
      },
    })] : []),
    ...(process.env.NODE_ENV === 'development' ? [anonymous()] : []),
    // API keys (created on the account page) act as their user, given in the `x-api-key` header
    //  or as `Authorization: Bearer <key>` (docs/user/api.md)
    apiKey({
      schema: { apikey: { modelName: 'auth_api_key' } },
      enableSessionForAPIKeys: true,
      defaultPrefix: 'pwb_',
      customAPIKeyGetter: (ctx) => {
        const header = ctx.headers?.get('x-api-key')
        if (header) return header
        const bearer = /^Bearer\s+(pwb_\S+)$/.exec(ctx.headers?.get('authorization') ?? '')
        return bearer ? bearer[1] : null
      },
      // (the default allows 10 requests a day, an upload alone takes several)
      rateLimit: { enabled: false },
      keyExpiration: { defaultExpiresIn: null, maxExpiresIn: 365 },
    }),
  ],
  advanced: {
    database: {
      generateId: 'uuid',
      validateSchema: false,
    },
  },
})

export type SessionWithId = {
  user: { id: string, name?: string | null, email?: string | null, image?: string | null }
  expires: string
}

