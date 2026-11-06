import jwtDecode from 'jwt-decode'
import type { OAuthConfig, OAuthUserConfig } from "next-auth/providers/index"
import type { DB } from '@/db'

const tokenColumns = ['access_token', 'refresh_token', 'expires_at', 'token_type', 'scope', 'id_token', 'session_state'] as const

export default function ORCIDProvider<P extends { sub: string, given_name: string, family_name: string, }>(options: OAuthUserConfig<P>, db: DB): OAuthConfig<P> {
  return {
    id: 'orcid',
    name: 'ORCID',
    type: 'oauth',
    wellKnown: "https://orcid.org/.well-known/openid-configuration",
    authorization: { params: { scope: "openid" } },
    idToken: true,
    checks: ["pkce", "state"],
    async profile(profile, tokens) {
      // keep account up to date
      if (tokens.id_token) {
        Object.assign(tokens, { expires_at: jwtDecode<{exp: string}>(tokens.id_token).exp })
      }
      const data = Object.fromEntries(tokenColumns.flatMap(key => key in tokens ? [[key, (tokens as Record<string, unknown>)[key] ?? null]] : []))
      if (Object.keys(data).length > 0) {
        await db.updateTable('account')
          .set(data)
          .where('provider', '=', 'orcid')
          .where('providerAccountId', '=', profile.sub)
          .execute()
      }
      return {
        id: profile.sub,
        name: `${profile.given_name} ${profile.family_name}`,
      }
    },
    style: {
      logo: 'https://orcid.org/assets/vectors/orcid.logo.icon.svg',
      logoDark: 'https://orcid.org/assets/vectors/orcid.logo.icon.svg',
      bgDark: "#fff",
      bg: "#fff",
      text: "#000",
      textDark: "#000",
    },
    options,
  }
}
