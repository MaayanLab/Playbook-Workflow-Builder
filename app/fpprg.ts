import krg from '@/app/krg'
import FPPRG from '@/core/FPPRG'
import db from '@/app/db'
import { start_workers } from '@/core/engine'
import cache from '@/utils/global_cache'

export default cache('fpprg', () => {
  const fpprg = new FPPRG(db)
  const n_workers = +(process.env.N_WORKERS||'1')
  if (n_workers > 0) {
    start_workers(krg, fpprg, n_workers)
  }
  return fpprg
})
