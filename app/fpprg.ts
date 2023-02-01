import krg from '@/app/krg'
import FPPRG from '@/core/FPPRG'
import db from '@/app/db'
import { start_workers } from '@/core/engine'
import cache from '@/utils/global_cache'

export default cache('fpprg', () => {
  const fpprg = new FPPRG(db)
  start_workers(krg, fpprg, +(process.env.N_WORKERS||'1'))
  return fpprg
})
