import krg from '@/app/krg'
import { PgDatabase, MemoryDatabase, Database } from '@/core/FPPRG'
import { process_insertion_dispatch } from '@/core/engine'
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
const fpprg = global.fpprg || (process.env.DATABASE_URL ? new PgDatabase(process.env.DATABASE_URL) : new MemoryDatabase())
if (global.fpprg && global.detach) global.detach()
const detach = process_insertion_dispatch(krg, fpprg)
global.fpprg = fpprg
global.detach = detach
export default fpprg
