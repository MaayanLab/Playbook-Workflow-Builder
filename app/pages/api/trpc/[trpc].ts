import { createNextApiHandler } from '@trpc/server/adapters/next'
import { appRouter } from '@/app/trpc/router'
import { createNextContext } from '@/app/trpc/context'

export default createNextApiHandler({
  router: appRouter,
  createContext: createNextContext,
  onError: ({ error, path }) => {
    if (error.code === 'INTERNAL_SERVER_ERROR') console.error(`[trpc] ${path}:`, error)
  },
})
