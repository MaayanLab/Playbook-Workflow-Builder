import handler from '@/utils/next-rest'
import db from '@/app/db'
import { NotFoundError } from '@/spec/error'
import { z } from 'zod'
import { fileAsStream } from "@/components/core/file/api/download"

const QueryType = z.object({
  object_id: z.string(),
  access_id: z.string(),
})
const BodyType = z.object({
  pasports: z.array(z.string()).optional(),
})

export default handler(async (req, res) => {
  const { object_id, access_id } = QueryType.parse(req.query)
  if (req.method === 'POST') {
    const body = BodyType.parse(req.body)
    // TODO: potentially authorize with passports
  }
  if (access_id !== 'https') throw new NotFoundError()
  const upload = await db.objects.user_upload_complete.findUnique({ where: { id: object_id } })
  if (upload === null) throw new NotFoundError()
  const fileStream = await fileAsStream(upload)
  fileStream.pipe(res)
})
