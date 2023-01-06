import PgBoss from 'pg-boss'
import KRG from '@/core/KRG'
import { PgDatabase } from '@/core/FPPRG'
import { resolve_process } from '@/core/engine'

/**
 * This worker receives jobs from the boss work-queue which ensures
 *  submitted jobs are given to one and only one worker. This should
 *  be run with several replicas in production.
 */
async function worker() {
  if (!process.env.DATABASE_URL) throw new Error('Missing `DATABASE_URL`')
  const boss = new PgBoss(process.env.DATABASE_URL)
  const db = new PgDatabase(process.env.DATABASE_URL)
  const krg = new KRG()
  boss.on('error', error => console.error(error))
  await boss.start()
  await boss.work('work-queue', async (job) => {
    // the job.data should contain the process id
    const processId = job.data as string
    console.log(`Processing ${processId}...`)
    // we fetch it from the db
    const instanceProcess = await db.getProcess(job.data as string)
    if (!instanceProcess) throw new Error(`Process ${job.data} not found`)
    if (instanceProcess.resolved === undefined) {
      // resolve the process
      const resolved = await resolve_process(krg, instanceProcess)
      // store the result in the db
      await db.upsertResolved(resolved)
    }
  })
}

worker().catch(error => console.error(error))
