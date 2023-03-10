import type { OAuthConfig, OAuthUserConfig } from "next-auth/providers"

export default function ORCIDProvider<P extends { sub: string, given_name: string, family_name: string, }>(options: OAuthUserConfig<P>): OAuthConfig<P> {
  return {
    id: 'orcid',
    name: 'ORCID',
    type: 'oauth',
    wellKnown: "https://orcid.org/.well-known/openid-configuration",
    authorization: { params: { scope: "openid" } },
    idToken: true,
    checks: ["pkce", "state"],
    profile(profile) {
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
