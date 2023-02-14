import type { OAuthConfig, OAuthUserConfig } from "next-auth/providers"

export default function ORCIDProvider<P extends { sub:string, name: string, email: string }>(options: OAuthUserConfig<P>): OAuthConfig<P> {
  return {
    id: 'orcid',
    name: 'ORCID',
    type: 'oauth',
    wellKnown: "https://orcid.org/.well-known/openid-configuration",
    authorization: { params: { scope: "openid email" } },
    idToken: true,
    checks: ["pkce", "state"],
    profile(profile) {
      return {
        id: profile.sub,
        name: profile.name,
        email: profile.email,
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
