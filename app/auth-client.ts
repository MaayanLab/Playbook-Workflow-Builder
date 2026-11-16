import { createAuthClient } from 'better-auth/react'
import { anonymousClient, magicLinkClient } from 'better-auth/client/plugins'
import { apiKeyClient } from '@better-auth/api-key/client'

export const authClient = createAuthClient({
  plugins: [anonymousClient(), magicLinkClient(), apiKeyClient()],
})

