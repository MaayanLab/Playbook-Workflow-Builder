import { createAuthClient } from 'better-auth/react'
import { anonymousClient, magicLinkClient } from 'better-auth/client/plugins'

export const authClient = createAuthClient({
  plugins: [anonymousClient(), magicLinkClient()],
})

