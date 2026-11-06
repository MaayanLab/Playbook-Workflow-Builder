import fpprg from '@/app/fpprg'
import { bus } from '@/app/db'
import { adminProcedure, router } from '../init'

export const adminRouter = router({
  jobs: adminProcedure
    .query(async () => {
      const jobs = await bus.jobs('work-queue')
      return await Promise.all(jobs.map(async (job) => {
        const proc = await fpprg.getProcess(job.id)
        if (!proc) return null
        return { ...job, type: proc.type }
      }))
    }),
})
