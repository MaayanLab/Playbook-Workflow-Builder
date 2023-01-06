import PgBoss from 'pg-boss'
import { PgDatabase, Process } from '@/core/FPPRG'

/**
 * This dispatcher submits jobs to the boss work-queue
 *  making it to the workers. Only one of these should be run.
 */
async function dispatcher() {
  if (!process.env.DATABASE_URL) throw new Error('Missing `DATABASE_URL`')
  const boss = new PgBoss(process.env.DATABASE_URL)
  const db = new PgDatabase(process.env.DATABASE_URL)
  boss.on('error', error => console.error(error))
  await boss.start()
  return db.listen(async (table, record) => {
    if (table === 'process') {
      const instanceProcess = record as Process
      if (instanceProcess.resolved === undefined) {
        await boss.send('work-queue', { id: instanceProcess.id })
      }
    }
  })
}

dispatcher().catch(error => console.error(error))
