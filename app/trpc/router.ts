import { createCallerFactory, router } from './init'
import { adminRouter } from './routers/admin'
import { chatRouter } from './routers/chat'
import { fpprgRouter } from './routers/fpprg'
import { integrationsRouter } from './routers/integrations'
import { playbooksRouter } from './routers/playbooks'
import { usersRouter } from './routers/users'
import { componentRouters } from '@/components/trpc'
import { systemRouter } from './routers/system'
import { catalogRouter } from './routers/catalog'
import { publicApiRouter } from './routers/public-api'
import { tunnelRouter } from './routers/tunnel'
import { exportsRouter } from './routers/exports'

export const appRouter = router({
  admin: adminRouter,
  chat: chatRouter,
  fpprg: fpprgRouter,
  integrations: integrationsRouter,
  playbooks: playbooksRouter,
  users: usersRouter,
  components: componentRouters,
  system: systemRouter,
  catalog: catalogRouter,
  publicApi: publicApiRouter,
  tunnel: tunnelRouter,
  exports: exportsRouter,
})

export type AppRouter = typeof appRouter

export const createCaller = createCallerFactory(appRouter)
