import krg from '@/app/krg'
import FPPRG from '@/core/FPPRG'
import db from '@/app/db'
import { start_workers } from '@/core/engine'
declare global {
  var fpprg: FPPRG | undefined
}

/**
 * This trick persists the database in a global object across hot-reloads
 * 
 * If the global objects exist, we'll detach the old engine
 *  and attach the new one to the persistent fpprg
 */
let fpprg: FPPRG
if (!global.fpprg) {
  global.fpprg = fpprg = new FPPRG(db)
  start_workers(krg, fpprg, +(process.env.N_WORKERS||'1'))
} else {
  fpprg = global.fpprg
}

export default fpprg
