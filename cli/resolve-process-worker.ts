import db from '@/app/db'
import fpprg from '@/app/fpprg'
import krg from '@/app/krg'
import { start_workers } from '@/core/engine'
import { start_semantic_search_worker } from '@/core/semantic-search/worker'

if (process.env.N_WORKERS && +process.env.N_WORKERS) {
  start_workers(krg, fpprg, +process.env.N_WORKERS)
  start_semantic_search_worker(db)
} else {
  throw new Error("N_WORKERS environment variable is required");
}
