import krg from '@/app/krg'
import { MemoryDatabase, Database } from '@/core/FPPRG'
import attach from '@/core/engine'
declare global {
  var fpprg: Database | undefined
  var detach: (() => void) | undefined
}

/**
 * This trick persists the in-memory database in
 *  a global object across hot-reloads
 * 
 * If the global objects exist, we'll detach the old engine
 *  and attach the new one to the persistent fpprg
 */
const fpprg = global.fpprg || new MemoryDatabase()
if (global.fpprg && global.detach) global.detach()
const detach = attach(krg, fpprg)
global.fpprg = fpprg
global.detach = detach
export default fpprg
