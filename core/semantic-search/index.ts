import type DB from "@/app/db"
import { NotFoundError } from "@/spec/error"
import { fpl_expand } from "@/core/common"
import FPPRG from "@/core/FPPRG"
import KRG from "@/core/KRG"
import OpenAI from "openai"

export async function indexSoon(props: { db: typeof DB, fpl_id: string }) {
    await props.db.send('embed', { id: props.fpl_id })
}
export async function index(props: { openai: OpenAI, db: typeof DB, fpprg: FPPRG, krg: KRG, fpl_id: string }) {
    const fpl = await props.fpprg.getFPL(props.fpl_id)
    if (!fpl) throw new NotFoundError()
    const fpl_embedding = await props.db.objects.fpl_embedding.findUnique({ where: { id: props.fpl_id } })
    if (fpl_embedding) return fpl_embedding.embedding
    const { story } = await fpl_expand({ krg: props.krg, fpl })
    const input = story.ast.flatMap(part => part.tags.includes('abstract') ? [part.type === 'bibitem' ? '\n' : '', part.text] : []).join('')
    const embeddings = await props.openai.embeddings.create({
        input,
        model: 'text-embedding-3-small',
        dimensions: 1536,
    })
    await props.db.objects.fpl_embedding.create({
        data: {
            id: props.fpl_id,
            embedding: embeddings.data[0].embedding,
        }
    })
    return embeddings.data[0].embedding
}
