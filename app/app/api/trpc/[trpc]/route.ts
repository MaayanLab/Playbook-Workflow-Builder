import { fetchRequestHandler } from '@trpc/server/adapters/fetch'
import { appRouter } from '@/app/trpc/router'
import { createFetchContext } from '@/app/trpc/context'

function handler(req: Request) {
  return fetchRequestHandler({
    endpoint: '/api/trpc',
    req,
    router: appRouter,
    createContext: () => createFetchContext(req),
    onError: ({ error, path }) => {
      if (error.code === 'INTERNAL_SERVER_ERROR') console.error(`[trpc] ${path}:`, error)
    },
  })
}

export { handler as GET, handler as POST }
