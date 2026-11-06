import { publicProcedure, router } from '../init'

export const systemRouter = router({
  runtimeConfig: publicProcedure
    .query(() => Object.fromEntries(Object.entries(process.env).filter(([key]) => key.startsWith('NEXT_PUBLIC_'))) as Record<string, string>),
})
