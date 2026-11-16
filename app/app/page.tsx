import { redirect } from 'next/navigation'

// the landing page is configured at runtime (LANDING_PAGE), so this can't be a redirect in next.config.js
export const dynamic = 'force-dynamic'

export default function Page() {
  const landingPage = process.env.LANDING_PAGE
  redirect(landingPage && landingPage !== '/' ? landingPage : '/graph/extend')
}
