import { Database } from '@/core/FPPRG'

declare global {
  var fpprg: Database | undefined
}

/**
 * This trick persists the in-memory database in
 *  a global object across hot-reloads
 */
const fpprg = global.fpprg || new Database()

if (process.env.NODE_ENV !== 'production') {
  global.fpprg = fpprg
}

export default fpprg
