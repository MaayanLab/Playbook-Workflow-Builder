import krg from '@/app/krg'
import { PgDatabase, MemoryDatabase, Database } from '@/core/FPPRG'
import { process_insertion_dispatch, start_workers } from '@/core/engine'
declare global {
  var fpprg: Database | undefined
  var detach: (() => void) | undefined
}

/**
 * This trick persists the database in a global object across hot-reloads
 * 
 * If the global objects exist, we'll detach the old engine
 *  and attach the new one to the persistent fpprg
 */
let fpprg: Database
if (process.env.DATABASE_URL) {
  if (!global.fpprg) {
    global.fpprg = fpprg = new PgDatabase(process.env.DATABASE_URL)
  } else {
    fpprg = global.fpprg
  }
  if (process.env.N_WORKERS && +process.env.N_WORKERS) {
    global.detach = start_workers(krg, fpprg as PgDatabase, +process.env.N_WORKERS)
  }
} else {
  if (!global.fpprg) {
    global.fpprg = fpprg = new MemoryDatabase()
  } else {
    fpprg = global.fpprg
  }
  if (global.detach) global.detach()
  global.detach = process_insertion_dispatch(krg, fpprg)
}

export default fpprg
