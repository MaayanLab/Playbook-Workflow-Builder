import type DB from "@/app/db";
import { z } from 'zod'
import openai from '@/app/extensions/openai'
import { fpl_expand } from "./common";
import krg from '@/app/krg'
import fpprg from '@/app/fpprg'
import { NotFoundError } from "@/spec/error";

const JobC = z.object({ data: z.object({ id: z.string() }) })

export function start_semantic_search_worker(db: typeof DB) {
    db.work('embed', { teamSize: 1, teamConcurrency: 1 }, async (job) => {
        const { data: { id } } = JobC.parse(job)
        const fpl = await fpprg.getFPL(id)
        if (!fpl) throw new NotFoundError()
        const fpl_embedding = await db.objects.fpl_embedding.findUnique({ where: { id } })
        if (fpl_embedding) return
        const { story } = await fpl_expand({ krg, fpl })
        const input = story.ast.flatMap(part => part.tags.includes('abstract') ? [part.type === 'bibitem' ? '\n' : '', part.text] : []).join('')
        const embeddings = await (await openai).embeddings.create({
            input,
            model: 'text-embedding-3-small',
            dimensions: 1536,
        })
        await db.objects.fpl_embedding.create({
            data: {
                id,
                embedding: embeddings.data[0].embedding,
            }
        })
    })
}
