import { getKRG } from '@/app/server-krg'
import FPPRG from '@/core/FPPRG'
import db, { bus } from '@/app/db'
import { start_workers } from '@/core/engine'
import { start_export_workers } from '@/app/exports'
import { start_dedupe_workers } from '@/components/core/file/server/dedupe'
import cache from '@/utils/global_cache'

export default cache('fpprg', () => {
  const fpprg = new FPPRG(db, bus)
  const n_workers = +(process.env.N_WORKERS||'1')
  if (n_workers > 0) {
    start_workers(getKRG, fpprg, n_workers)
    // exports (LaTeX, PDFs, ...) are produced by the same workers, in a separate queue so
    //  that long running exports don't hold up process resolution
    start_export_workers(getKRG, fpprg, +(process.env.N_EXPORT_WORKERS||'1'))
    // & browser uploads are verified & deduplicated (components/core/file/server/dedupe.ts)
    start_dedupe_workers(1)
  }
  return fpprg
})
