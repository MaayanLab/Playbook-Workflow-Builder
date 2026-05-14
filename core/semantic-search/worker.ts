import type DB from "@/app/db";
import { z } from 'zod'
import openai from '@/app/extensions/openai'
import krg from '@/app/krg'
import fpprg from '@/app/fpprg'
import { index } from ".";

const JobC = z.object({ data: z.object({ id: z.string() }) })

export function start_semantic_search_worker(db: typeof DB) {
  db.work('embed', { teamSize: 1, teamConcurrency: 1 }, async (job) => {
    const { data: { id } } = JobC.parse(job)
    await index({
        openai: await openai,
        db,
        fpprg,
        krg,
        fpl_id: id,
    })
  })
}
