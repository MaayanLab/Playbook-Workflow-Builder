import { createCallerFactory, router } from './init'
import { adminRouter } from './routers/admin'
import { chatRouter } from './routers/chat'
import { fpprgRouter } from './routers/fpprg'
import { integrationsRouter } from './routers/integrations'
import { playbooksRouter } from './routers/playbooks'

export const appRouter = router({
  admin: adminRouter,
  chat: chatRouter,
  fpprg: fpprgRouter,
  integrations: integrationsRouter,
  playbooks: playbooksRouter,
})

export type AppRouter = typeof appRouter

export const createCaller = createCallerFactory(appRouter)
