import fpprg from '@/app/fpprg'
import krg from '@/app/krg'
import { start_workers } from '@/core/engine'

if (process.env.N_WORKERS && +process.env.N_WORKERS) {
  start_workers(krg, fpprg, +process.env.N_WORKERS)
} else {
  throw new Error("N_WORKERS environment variable is required");
}
