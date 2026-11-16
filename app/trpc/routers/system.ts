import { publicProcedure, router } from '../init'
import { getRuntimeConfig } from '@/app/runtime-config'

export const systemRouter = router({
  runtimeConfig: publicProcedure
    .query(() => getRuntimeConfig()),
})
