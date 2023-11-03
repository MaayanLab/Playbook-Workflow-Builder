import handler from '@/utils/next-rest'
import db from '@/app/db'
import { NotFoundError } from '@/spec/error'
import { z } from 'zod'

const QueryType = z.object({
  object_id: z.string(),
})
const BodyType = z.object({
  expand: z.boolean().optional(),
  passports: z.array(z.string()).optional(),
})

export default handler(async (req, res) => {
  const { object_id } = QueryType.parse(req.query)
  if (req.method === 'POST') {
    const body = BodyType.parse(req.body)
    // TODO: potentially authorize with passports
  }
  const upload = await db.objects.user_upload_complete.findUnique({ where: { id: object_id } })
  if (upload === null) throw new NotFoundError()
  if (req.method === 'OPTIONS') {
    res.status(200)
    return
  }
  res.json({
    "id": upload.id,
    "name": upload.filename,
    "self_uri": `${(process.env.PUBLIC_URL||'').replace(/^https?:/, 'drs:').replace(/\/$/, '')}/${upload.id}`,
    "size": upload.size,
    "created_time": upload.created,
    "checksums": [
      {"type": "sha-256", "checksum": upload.sha256},
    ],
    "access_methods": [
      {'type': 'https', 'access_id': 'https'},
      {'type': 'https', 'access_url': `${process.env.PUBLIC_URL||''}/ga4gh/drs/v1/objects/${object_id}/access/https/data`},
    ],
  })
})
